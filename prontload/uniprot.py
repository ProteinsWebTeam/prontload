#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import math
import os
import tempfile


def get_comments(cur):
    cur.execute(
        """
        SELECT E.ACCESSION, B.COMMENT_TOPICS_ID, T.TOPIC, SS.TEXT
        FROM SPTR.DBENTRY@SWPREAD E
        INNER JOIN SPTR.COMMENT_BLOCK@SWPREAD B ON E.DBENTRY_ID = B.DBENTRY_ID
        INNER JOIN SPTR.CV_COMMENT_TOPICS@SWPREAD T ON B.COMMENT_TOPICS_ID = T.COMMENT_TOPICS_ID
        INNER JOIN SPTR.COMMENT_STRUCTURE@SWPREAD S ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
        INNER JOIN SPTR.COMMENT_SUBSTRUCTURE@SWPREAD SS ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
        WHERE E.ENTRY_TYPE = 0
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        """
    )

    topics = {}
    protein_comments = set()
    for row in cur:
        protein_ac = row[0]
        topic_id = math.floor(row[1])
        comment = row[3]

        if topic_id not in topics:
            topics[topic_id] = {
                'name': row[2],
                'comments': {}
            }

        comments = topics[topic_id]['comments']
        if comment in comments:
            comment_id = comments[comment]
        else:
            comment_id = comments[comment] = len(comments) + 1

        protein_comments.add((protein_ac, topic_id, comment_id))

    comments = []
    for topic_id, topic in topics.items():
        for comment, comment_id in topic['comments'].items():
            comments.append((topic_id, comment_id, comment))

    topics = [(topic_id, topic['name']) for topic_id, topic in topics.items()]

    return topics, comments, list(protein_comments)


def _dump_descriptions(descriptions, tmpdir=None):
    fd, filepath = tempfile.mkstemp(suffix='.json.gz', dir=tmpdir)
    os.close(fd)

    with gzip.open(filepath, 'wt') as fh:
        json.dump([{'text': descr, 'proteins': descriptions[descr]} for descr in sorted(descriptions)], fh)

    return filepath


def dump_descriptions(cur, tmpdir=None, chunk_size=1000000):
    if tmpdir and not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    cur.execute(
        """
        SELECT DESCR, ACCESSION
        FROM (
          SELECT
            E.ACCESSION,
            D.DESCR,
            ROW_NUMBER() OVER (PARTITION BY E.ACCESSION ORDER BY D.SECTION_GROUP_ID, D.DESC_ID) R
          FROM SPTR.DBENTRY@SWPREAD E
            LEFT OUTER JOIN SPTR.DBENTRY_2_DESC@SWPREAD D ON E.DBENTRY_ID = D.DBENTRY_ID
          WHERE E.ENTRY_TYPE IN (0, 1)
                AND E.MERGE_STATUS != 'R'
                AND E.DELETED = 'N'
                AND E.FIRST_PUBLIC IS NOT NULL
        )
        WHERE R = 1
        ORDER BY DESCR
        """
    )

    descriptions = {}
    files = []
    for descr, accession in cur:
        if descr in descriptions:
            d = descriptions[descr]
        else:
            if len(descriptions) == chunk_size:
                files.append(_dump_descriptions(descriptions, tmpdir=tmpdir))
                descriptions = {}

            d = descriptions[descr] = []

        d.append(accession)

    if descriptions:
        files.append(_dump_descriptions(descriptions, tmpdir=tmpdir))

    return files
