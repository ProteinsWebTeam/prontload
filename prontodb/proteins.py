#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import logging
import math
import os
import tempfile

from . import oracledb


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


def get_descriptions(cur, tmpdir=None, chunk_size=1000000):
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


def load_comments(dsn, schema, **kwargs):
    chunk_size = kwargs.get('chunk_size', 100000)

    logging.info('loading UniProt comments')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'CV_COMMENT_TOPIC')
    oracledb.drop_table(cur, schema, 'COMMENT_VALUE')
    oracledb.drop_table(cur, schema, 'PROTEIN_COMMENT')

    cur.execute(
        """
        CREATE TABLE {}.CV_COMMENT_TOPIC
        (
            TOPIC_ID NUMBER(2) NOT NULL,
            TOPIC VARCHAR2(30) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        CREATE TABLE {}.COMMENT_VALUE
        (
            TOPIC_ID NUMBER(2) NOT NULL,
            COMMENT_ID NUMBER(6) NOT NULL,
            TEXT VARCHAR2(4000) NOT NULL,
            CONSTRAINT PK_COMMENT_VALUE PRIMARY KEY (TOPIC_ID, COMMENT_ID)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        CREATE TABLE {}.PROTEIN_COMMENT
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            TOPIC_ID NUMBER(2) NOT NULL,
            COMMENT_ID NUMBER(6) NOT NULL,
            CONSTRAINT PK_PROTEIN_COMMENT PRIMARY KEY (PROTEIN_AC, TOPIC_ID, COMMENT_ID)
        ) NOLOGGING
        """.format(schema)
    )

    topics, comments, protein_comments = get_comments(cur)

    cur.executemany(
        """
        INSERT /*+APPEND*/ INTO {}.CV_COMMENT_TOPIC (TOPIC_ID, TOPIC)
        VALUES (:1, :2)
        """.format(schema),
        topics
    )
    con.commit()

    for i in range(0, len(comments), chunk_size):
        cur.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.COMMENT_VALUE (TOPIC_ID, COMMENT_ID, TEXT)
            VALUES (:1, :2, :3)
            """.format(schema),
            comments[i:i + chunk_size]
        )
        con.commit()

    for i in range(0, len(protein_comments), chunk_size):
        cur.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.PROTEIN_COMMENT (PROTEIN_AC, TOPIC_ID, COMMENT_ID)
            VALUES (:1, :2, :3)
            """.format(schema),
            protein_comments[i:i + chunk_size]
        )
        con.commit()

    oracledb.gather_stats(cur, schema, 'CV_COMMENT_TOPIC', cascade=True)
    oracledb.gather_stats(cur, schema, 'COMMENT_VALUE', cascade=True)
    oracledb.gather_stats(cur, schema, 'PROTEIN_COMMENT', cascade=True)

    oracledb.grant(cur, 'SELECT', schema, 'CV_COMMENT_TOPIC', 'INTERPRO_SELECT')
    oracledb.grant(cur, 'SELECT', schema, 'COMMENT_VALUE', 'INTERPRO_SELECT')
    oracledb.grant(cur, 'SELECT', schema, 'PROTEIN_COMMENT', 'INTERPRO_SELECT')

    cur.close()
    con.close()
    logging.info('\tcomplete')


def load_descriptions(dsn, schema, **kwargs):
    chunk_size = kwargs.get('chunk_size', 100000)
    tmpdir = kwargs.get('tmpdir', tempfile.gettempdir())

    logging.info('loading UniProt descriptions')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'DESC_VALUE')
    oracledb.drop_table(cur, schema, 'PROTEIN_DESC')

    cur.execute(
        """
        CREATE TABLE {}.DESC_VALUE
        (
            DESC_ID NUMBER(10) NOT NULL,
            TEXT VARCHAR2(4000) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        CREATE TABLE {}.PROTEIN_DESC
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            DESC_ID NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    files = get_descriptions(cur, chunk_size=chunk_size, tmpdir=tmpdir)

    desc_id = 1
    for i, filepath in enumerate(files):

        with gzip.open(filepath, 'rt') as fh:
            descriptions = json.load(fh)

        os.unlink(filepath)

        cv_table = []
        rel_table = []
        for descr in descriptions:
            cv_table.append((desc_id, descr['text']))
            rel_table += [(accession, desc_id) for accession in descr['proteins']]
            desc_id += 1

        cur.executemany(
            """
            INSERT /*+ APPEND */ INTO {}.DESC_VALUE (DESC_ID, TEXT)
            VALUES (:1, :2)
            """.format(schema),
            cv_table
        )
        # Commit after each transaction to avoid ORA-12838
        con.commit()

        for j in range(0, len(rel_table), chunk_size):
            cur.executemany(
                """
                INSERT /*+ APPEND */ INTO {}.PROTEIN_DESC (PROTEIN_AC, DESC_ID)
                VALUES (:1, :2)
                """.format(schema),
                rel_table[j:j + chunk_size]
            )
            # Commit after each transaction to avoid ORA-12838
            con.commit()

    cur.execute(
        """
        ALTER TABLE {}.DESC_VALUE
        ADD CONSTRAINT PK_DESC_VALUE PRIMARY KEY (DESC_ID)
        """.format(schema)
    )

    cur.execute(
        """
        ALTER TABLE {}.PROTEIN_DESC
        ADD CONSTRAINT PK_PROTEIN_DESC PRIMARY KEY (PROTEIN_AC, DESC_ID)
        """.format(schema)
    )

    oracledb.gather_stats(cur, schema, 'DESC_VALUE', cascade=True)
    oracledb.gather_stats(cur, schema, 'PROTEIN_DESC', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'DESC_VALUE', 'INTERPRO_SELECT')
    oracledb.grant(cur, 'SELECT', schema, 'PROTEIN_DESC', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('\tcomplete')


def load_enzymes(dsn, schema, **kwargs):
    logging.info('loading ENZYME annotations')

    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'ENZYME')

    cur.execute(
        """
        CREATE TABLE {}.ENZYME
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ECNO VARCHAR2(15) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.ENZYME (PROTEIN_AC, ECNO)
        SELECT DISTINCT E.ACCESSION, D.DESCR
        FROM SPTR.DBENTRY@SWPREAD E
        LEFT OUTER JOIN SPTR.DBENTRY_2_DESC@SWPREAD D ON E.DBENTRY_ID = D.DBENTRY_ID
        LEFT OUTER JOIN SPTR.CV_DESC@SWPREAD C ON D.DESC_ID = C.DESC_ID
        WHERE E.ENTRY_TYPE IN (0, 1)
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        AND C.SUBCATG_TYPE='EC'
        """.format(schema)
    )
    con.commit()

    cur.execute('CREATE INDEX I_ENZYME$PROTEIN ON {}.ENZYME (PROTEIN_AC) NOLOGGING'.format(schema))
    cur.execute('CREATE INDEX I_ENZYME$EC ON {}.ENZYME (ECNO) NOLOGGING'.format(schema))
    oracledb.gather_stats(cur, schema, 'ENZYME', cascade=True)

    oracledb.grant(cur, 'SELECT', schema, 'ENZYME', 'INTERPRO_SELECT')

    cur.close()
    con.close()

    logging.info('\tcomplete')


def load_proteins(dsn, schema, **kwargs):
    logging.info('loading proteins')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'PROTEIN')

    cur.execute(
        """
        CREATE TABLE {}.PROTEIN
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            NAME VARCHAR2(16) NOT NULL ,
            DBCODE CHAR(1) NOT NULL,
            LEN NUMBER(5) NOT NULL,
            FRAGMENT CHAR(1) NOT NULL,
            TAX_ID NUMBER(15)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.PROTEIN (PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, TAX_ID)
        SELECT PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, TAX_ID
        FROM INTERPRO.PROTEIN
        """.format(schema)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.PROTEIN
        ADD CONSTRAINT PK_PROTEIN PRIMARY KEY (PROTEIN_AC)
        """.format(schema)
    )
    cur.execute('CREATE INDEX I_PROTEIN$DBCODE ON {}.PROTEIN (DBCODE) NOLOGGING'.format(schema))
    cur.execute('CREATE INDEX I_PROTEIN$NAME ON {}.PROTEIN (NAME) NOLOGGING'.format(schema))
    oracledb.gather_stats(cur, schema, 'PROTEIN', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'PROTEIN', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('\tcomplete')
