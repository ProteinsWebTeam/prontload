#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

from . import oracledb


def drop_all(dsn, schema, forgive_busy=False):
    con = oracledb.connect(dsn)
    cur = con.cursor()
    tables = oracledb.list_tables(cur, schema)
    for table in tables:
        if oracledb.drop_table(cur, schema, table, forgive_busy=forgive_busy):
            logging.info('table {} dropped'.format(table))
        else:
            logging.error('table {} could not be dropped'.format(table))

    cur.close()
    con.close()


def create_synonyms(dsn, schema, **kwargs):
    logging.info('creating synonyms')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    for obj in ('ENTRY', 'ENTRY2METHOD', 'ENTRY2ENTRY', 'ENTRY2COMP'):
        oracledb.create_synonym(cur, 'INTERPRO', schema, obj)
    cur.close()
    con.close()
    logging.info('creating synonyms: complete')


def load_databases(dsn, schema, **kwargs):
    logging.info('loading databases')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'CV_DATABASE')

    cur.execute(
        """
        CREATE TABLE {}.CV_DATABASE
        (
            DBCODE VARCHAR2(10) NOT NULL,
            DBNAME VARCHAR2(50) NOT NULL,
            DBSHORT VARCHAR2(10) NOT NULL,
            VERSION VARCHAR2(20),
            FILE_DATE DATE,
            CONSTRAINT PK_DATABASE PRIMARY KEY (DBCODE)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.CV_DATABASE (DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE)
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        """.format(schema)
    )
    con.commit()

    oracledb.gather_stats(cur, schema, 'CV_DATABASE', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'CV_DATABASE', 'INTERPRO_SELECT')
    cur.close()
    con.close()
    logging.info('loading databases: complete')


def load_taxonomies(dsn, schema, **kwargs):
    chunk_size = kwargs.get('chunk_size', 100000)

    logging.info('loading taxons')

    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'ETAXI')

    cur.execute(
        """
        CREATE TABLE {}.ETAXI
        (
            TAX_ID NUMBER(10) NOT NULL,
            PARENT_ID NUMBER(10),
            SCIENTIFIC_NAME VARCHAR2(255) NOT NULL,
            RANK VARCHAR2(50),
            LEFT_NUMBER NUMBER,
            RIGHT_NUMBER NUMBER,
            FULL_NAME VARCHAR2(513)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.ETAXI (TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, LEFT_NUMBER, RIGHT_NUMBER, FULL_NAME)
        SELECT TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, LEFT_NUMBER, RIGHT_NUMBER, FULL_NAME
        FROM INTERPRO.ETAXI
        """.format(schema)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.ETAXI
        ADD CONSTRAINT PK_ETAXI PRIMARY KEY (TAX_ID)
        """.format(schema)
    )
    oracledb.gather_stats(cur, schema, 'ETAXI', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'ETAXI', 'INTERPRO_SELECT')
    logging.info('loading taxons: complete')

    logging.info('loading lineages')
    cur.execute('SELECT TAX_ID, PARENT_ID, LEFT_NUMBER, RANK FROM {}.ETAXI'.format(schema))
    taxons = {}
    for tax_id, parent_id, left_number, rank in cur:
        if tax_id == 131567:
            """
            taxID 131567 (cellular organisms) contains three superkingdoms:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)

            therefore it is not needed (we don't want a meta-superkingdom)
            """
            continue

        taxons[tax_id] = {
            'id': tax_id,
            'parent': parent_id,
            'left_number': left_number,
            'rank': 'superkingdom' if parent_id == 1 else rank
        }

    lineage = []
    for tax_id in taxons:

        t = taxons[tax_id]
        left_number = t['left_number']

        if not left_number:
            continue
        elif t['rank'] != 'no rank':
            lineage.append((left_number, t['id'], t['rank']))

        parent_id = t['parent']

        while parent_id:
            if parent_id not in taxons:
                # taxID 131567 missing from dictionary
                break

            t = taxons[parent_id]

            if t['rank'] != 'no rank':
                lineage.append((left_number, t['id'], t['rank']))

            parent_id = t['parent']

    oracledb.drop_table(cur, schema, 'LINEAGE')

    cur.execute(
        """
        CREATE TABLE {}.LINEAGE
        (
            LEFT_NUMBER NUMBER NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50)
        ) NOLOGGING
        """.format(schema)
    )

    for i in range(0, len(lineage), chunk_size):
        cur.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.LINEAGE (LEFT_NUMBER, TAX_ID, RANK)
            VALUES (:1, :2, :3)
            """.format(schema),
            lineage[i:i+chunk_size]
        )
        con.commit()

    cur.execute(
        """
        CREATE INDEX I_LINEAGE$L$R
        ON {}.LINEAGE (LEFT_NUMBER, RANK)
        NOLOGGING
        """.format(schema)
    )

    oracledb.gather_stats(cur, schema, 'LINEAGE', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'LINEAGE', 'INTERPRO_SELECT')

    cur.close()
    con.close()

    logging.info('loading lineages: complete')
