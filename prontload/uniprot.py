#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import os
import math
from tempfile import mkstemp

from .oracledb import Connection


def load_comments(dsn, schema, chunk_size=100000):
    tables = ("CV_COMMENT_TOPIC", "COMMENT_VALUE", "PROTEIN_COMMENT")
    con = Connection(dsn)
    for table in tables:
        con.drop_table(schema, table)

    con.execute(
        """
        CREATE TABLE {}.CV_COMMENT_TOPIC
        (
            TOPIC_ID NUMBER(2) NOT NULL,
            TOPIC VARCHAR2(30) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    con.execute(
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

    con.execute(
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

    query = """
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

    topics = {}
    protein2comment = set()
    for row in con.get(query):
        accession = row[0]
        topic_id = math.floor(row[1])
        comment = row[3]

        if topic_id not in topics:
            topics[topic_id] = {
                "name": row[2],
                "comments": {}
            }

        topic_comments = topics[topic_id]["comments"]
        if comment in topic_comments:
            comment_id = topic_comments[comment]
        else:
            comment_id = topic_comments[comment] = len(topic_comments) + 1

        protein2comment.add((accession, topic_id, comment_id))
    protein2comment = list(protein2comment)

    comments = []
    for topic_id, topic in topics.items():
        for comment, comment_id in topic["comments"].items():
            comments.append((topic_id, comment_id, comment))

    topics = [
        (topic_id, topic["name"])
        for topic_id, topic in topics.items()
    ]

    query = """
        INSERT /*+APPEND*/ INTO {}.CV_COMMENT_TOPIC (TOPIC_ID, TOPIC)
        VALUES (:1, :2)
    """.format(schema)
    con.executemany(query, topics)
    con.commit()

    query = """
        INSERT /*+APPEND*/ INTO {}.COMMENT_VALUE (TOPIC_ID, COMMENT_ID, TEXT)
        VALUES (:1, :2, :3)    
    """.format(schema)
    for i in range(0, len(comments), chunk_size):
        con.executemany(query, comments[i:i+chunk_size])
        con.commit()

    query = """
        INSERT /*+APPEND*/ INTO {}.PROTEIN_COMMENT (PROTEIN_AC, TOPIC_ID, COMMENT_ID)
        VALUES (:1, :2, :3)    
    """.format(schema)
    for i in range(0, len(protein2comment), chunk_size):
        con.executemany(query, protein2comment[i:i+chunk_size])
        con.commit()

    for table in tables:
        con.optimize_table(schema, table, cascade=True)
        con.grant("SELECT", schema, table, "INTERPRO_SELECT")


def load_descriptions(dsn, schema, tmpdir=None, chunk_size=100000):
    tables = ["DESC_VALUE", "PROTEIN_DESC"]
    con = Connection(dsn)
    for table in tables:
        con.drop_table(schema, table)

    con.execute(
        """
        CREATE TABLE {}.DESC_VALUE
        (
            DESC_ID NUMBER(10) NOT NULL,
            TEXT VARCHAR2(4000) NOT NULL
        ) NOLOGGING        
        """.format(schema)
    )

    con.execute(
        """
        CREATE TABLE {}.PROTEIN_DESC
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            DESC_ID NUMBER(10) NOT NULL
        ) NOLOGGING        
        """.format(schema)
    )

    query = """
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
    descriptions = []
    proteins = set()
    files = []
    _text = None
    for text, accession in con.get(query):
        if text != _text:
            if _text:
                descriptions.append({
                    "text": _text,
                    "proteins": list(proteins)
                })

                if len(descriptions) == chunk_size:
                    files.append(_dump_descriptions(descriptions, tmpdir))
                    descriptions = []

            proteins = set()
            _text = text

        proteins.add(accession)

    if proteins:
        descriptions.append({
            "text": _text,
            "proteins": list(proteins)
        })
        files.append(_dump_descriptions(descriptions, tmpdir))

    desc_id = 1
    for filepath in files:
        with gzip.open(filepath, "rt") as fh:
            descriptions = json.load(fh)

        os.unlink(filepath)

        cv_table = []
        rel_table = []
        for d in descriptions:
            cv_table.append((desc_id, d["text"]))
            rel_table += [(accession, desc_id) for accession in d["proteins"]]
            desc_id += 1

        con.executemany(
            """
            INSERT /*+ APPEND */ INTO {}.DESC_VALUE (DESC_ID, TEXT)
            VALUES (:1, :2)            
            """.format(schema),
            cv_table
        )
        # Commit after each transaction to avoid ORA-12838
        con.commit()

        for i in range(0, len(rel_table), chunk_size):
            con.executemany(
                """
                INSERT /*+ APPEND */ INTO {}.PROTEIN_DESC (PROTEIN_AC, DESC_ID)
                VALUES (:1, :2)                
                """.format(schema),
                rel_table[i:i+chunk_size]
            )
            con.commit()

    con.execute(
        """
        ALTER TABLE {}.DESC_VALUE
        ADD CONSTRAINT PK_DESC_VALUE PRIMARY KEY (DESC_ID)
        """.format(schema)
    )

    con.execute(
        """
        ALTER TABLE {}.PROTEIN_DESC
        ADD CONSTRAINT PK_PROTEIN_DESC PRIMARY KEY (PROTEIN_AC, DESC_ID)
        """.format(schema)
    )

    for table in tables:
        con.optimize_table(schema, table, cascade=True)
        con.grant("SELECT", schema, table, "INTERPRO_SELECT")


def _dump_descriptions(descriptions, tmpdir=None):
    fd, filepath = mkstemp(suffix=".json.gz", dir=tmpdir)
    os.close(fd)

    with gzip.open(filepath, "wt") as fh:
        json.dump(descriptions, fh)

    return filepath


def load_enzymes(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "ENZYME")
    con.execute(
        """
        CREATE TABLE {}.ENZYME
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            ECNO VARCHAR2(15) NOT NULL
        ) NOLOGGING        
        """.format(schema)
    )

    con.execute(
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

    con.execute(
        """
        CREATE INDEX I_ENZYME$PROTEIN 
        ON {}.ENZYME (PROTEIN_AC) NOLOGGING
        """.format(schema)
    )

    con.execute(
        """
        CREATE INDEX I_ENZYME$EC 
        ON {}.ENZYME (ECNO) NOLOGGING
        """.format(schema)
    )

    con.optimize_table(schema, "ENZYME", cascade=True)
    con.grant("SELECT", schema, "ENZYME", "INTERPRO_SELECT")