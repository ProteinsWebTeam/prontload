#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import logging
import math
import os
import tempfile
import time
from multiprocessing import Queue

from . import processor
from .. import oracledb, uniprot


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)


def drop_all(dsn, schema):
    con = oracledb.connect(dsn)
    cur = con.cursor()
    tables = oracledb.list_tables(cur, schema)
    for table in tables:
        logging.info('dropping {}'.format(table))
        oracledb.drop_table(cur, schema, table)

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
    logging.info('\tcomplete')


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

    topics, comments, protein_comments = uniprot.get_comments(cur)

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

    files = uniprot.dump_descriptions(cur, chunk_size=chunk_size, tmpdir=tmpdir)

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
    oracledb.gather_stats(cur, schema, 'PROTEIN', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'PROTEIN', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('\tcomplete')


def load_protein2go(dsn, schema, **kwargs):
    logging.info('loading protein GO annotations')

    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'PROTEIN2GO')

    cur.execute(
        """
        CREATE TABLE {}.PROTEIN2GO
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            GO_ID VARCHAR2(10) NOT NULL,
            EVIDENCE VARCHAR2(100) NOT NULL,
            REF_DB_CODE VARCHAR2(10),
            REF_DB_ID VARCHAR2(60)
        ) NOLOGGING
        """.format(schema)
    )

    """
    Filtering on length:
    Some annotations are not on proteins, but on post-translation modifications or processing events.
        e.g. P27958:PRO_0000037566 (protein: P27958; chain: PRO_0000037573)

    Protein accessions are 15 characters long (max), so anything above 15 characters cannot be an accession.
    A better (but heavier) approach would be to join with our PROTEIN table
    """
    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.PROTEIN2GO (PROTEIN_AC, GO_ID, EVIDENCE, REF_DB_CODE, REF_DB_ID)
        SELECT A.ENTITY_ID, A.GO_ID, E.GO_EVIDENCE, A.REF_DB_CODE, A.REF_DB_ID
        FROM GO.ANNOTATIONS@GOAPRO A
        INNER JOIN GO.ECO2EVIDENCE@GOAPRO E ON A.ECO_ID = E.ECO_ID
        INNER JOIN GO.CV_SOURCES@GOAPRO S ON S.CODE = A.SOURCE
        WHERE A.ENTITY_TYPE = 'protein'
        AND LENGTH(A.ENTITY_ID) <= 15
        AND E.GO_EVIDENCE != 'IEA'
        AND S.IS_PUBLIC = 'Y'
        ORDER BY A.ENTITY_ID
        """.format(schema)
    )
    con.commit()

    cur.execute('CREATE INDEX I_PROTEIN2GO$P$G ON {}.PROTEIN2GO (PROTEIN_AC, GO_ID) NOLOGGING'.format(schema))
    cur.execute('CREATE INDEX I_PROTEIN2GO$E ON {}.PROTEIN2GO (EVIDENCE) NOLOGGING'.format(schema))
    cur.execute('CREATE INDEX I_PROTEIN2GO$RC ON {}.PROTEIN2GO (REF_DB_CODE) NOLOGGING'.format(schema))
    oracledb.gather_stats(cur, schema, 'PROTEIN2GO', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'PROTEIN2GO', 'INTERPRO_SELECT')
    logging.info('\tcomplete')

    logging.info('loading publications')
    oracledb.drop_table(cur, schema, 'PUBLICATION')
    cur.execute(
        """
        CREATE TABLE {}.PUBLICATION
        (
            ID VARCHAR2(25) NOT NULL,
            TITLE VARCHAR2(1500),
            FIRST_PUBLISHED_DATE DATE
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.PUBLICATION (ID, TITLE, FIRST_PUBLISHED_DATE)
        SELECT ID, TITLE, FIRST_PUBLISH_DATE
        FROM GO.PUBLICATIONS@GOAPRO
        WHERE ID IN (
          SELECT DISTINCT A.REF_DB_ID
          FROM GO.ANNOTATIONS@GOAPRO A
            INNER JOIN GO.ECO2EVIDENCE@GOAPRO E ON A.ECO_ID = E.ECO_ID
            INNER JOIN GO.CV_SOURCES@GOAPRO S ON S.CODE = A.SOURCE
          WHERE A.ENTITY_TYPE = 'protein'
                AND LENGTH(A.ENTITY_ID) <= 15
                AND E.GO_EVIDENCE != 'IEA'
                AND S.IS_PUBLIC = 'Y'
                AND A.REF_DB_CODE = 'PMID'
        )
        """.format(schema)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.PUBLICATION
        ADD CONSTRAINT PK_PUBLICATION PRIMARY KEY (ID)
        """.format(schema)
    )
    oracledb.gather_stats(cur, schema, 'PUBLICATION', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'PUBLICATION', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('\tcomplete')


def load_methods(dsn, schema, **kwargs):
    logging.info('loading member database signatures')

    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'METHOD')

    cur.execute(
        """
        CREATE TABLE {}.METHOD
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            CANDIDATE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(220),
            SIG_TYPE CHAR(1)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.METHOD (METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE)
        SELECT METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE
        FROM INTERPRO.METHOD
        """.format(schema)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.METHOD
        ADD CONSTRAINT PK_METHOD PRIMARY KEY (METHOD_AC)
        """.format(schema)
    )
    cur.execute('CREATE INDEX I_METHOD$DBCODE ON {}.METHOD (DBCODE) NOLOGGING'.format(schema))
    oracledb.gather_stats(cur, schema, 'METHOD', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'METHOD', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('\tcomplete')


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
    logging.info('\tcomplete')

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

    logging.info('\tcomplete')


def load_terms(dsn, schema, **kwargs):
    logging.info('loading GO terms')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'TERM')

    cur.execute(
        """
        CREATE TABLE {}.TERM
        (
            GO_ID VARCHAR2(10) NOT NULL,
            NAME VARCHAR2(200) NOT NULL,
            CATEGORY CHAR(1) NOT NULL ,
            IS_OBSOLETE CHAR(1) NOT NULL ,
            DEFINITION VARCHAR2(4000),
            REPLACED_BY VARCHAR2(10)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.TERM (GO_ID, NAME, CATEGORY, IS_OBSOLETE, DEFINITION, REPLACED_BY)
        SELECT T.GO_ID GO_ID, T.NAME, T.CATEGORY, T.IS_OBSOLETE, D.DEFINITION, NULL
        FROM GO.TERMS@GOAPRO T
        INNER JOIN GO.DEFINITIONS@GOAPRO D ON T.GO_ID = D.GO_ID
        UNION ALL
        SELECT S.SECONDARY_ID GO_ID, T.NAME, T.CATEGORY, T.IS_OBSOLETE, D.DEFINITION, T.GO_ID
        FROM GO.SECONDARIES@GOAPRO S
        INNER JOIN GO.TERMS@GOAPRO T ON S.GO_ID = T.GO_ID
        INNER JOIN GO.DEFINITIONS@GOAPRO D ON T.GO_ID = D.GO_ID
        ORDER BY T.GO_ID
        """.format(schema)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.TERM
        ADD CONSTRAINT PK_TERM PRIMARY KEY (GO_ID)
        """.format(schema)
    )
    oracledb.gather_stats(cur, schema, 'TERM', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'TERM', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('\tcomplete')


def process_matches(dsn, schema, **kwargs):
    max_gap = kwargs.get('max_gap', 20)
    processes = kwargs.get('processes', 3)
    chunk_size = kwargs.get('chunk_size', 100000)
    max_size = kwargs.get('max_size', 10000000)

    logging.info('processing InterPro matches')

    proteins = Queue(maxsize=100000)
    comparisons = Queue(maxsize=100000)

    comparators = [
        processor.MatchComparator(proteins, comparisons, max_gap)
        for _ in range(min(1, processes-2))
    ]

    aggregator = processor.Aggregator(comparisons, dsn, schema, **kwargs)

    for p in comparators:
        p.start()

    aggregator.start()

    """
    MobiDB-lite: keep matches for MATCH table, but not for predictions

    PANTHER & PRINTS:
        Merge protein matches.
        If the signature is a family*, use the entire protein.

        * all PANTHER signatures, almost all PRINTS signatures
    """
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'MATCH')
    cur.execute(
        """
        CREATE TABLE {}.MATCH
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            STATUS CHAR(1) NOT NULL,
            DBCODE CHAR(1) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )
    cur.execute(
        """
        SELECT
          MA.PROTEIN_AC,
          MA.METHOD_AC,
          MA.POS_FROM,
          MA.POS_TO,
          MA.STATUS,
          MA.DBCODE,
          ME.SIG_TYPE,
          P.LEN,
          P.FRAGMENT
        FROM INTERPRO.MATCH MA
          INNER JOIN INTERPRO.METHOD ME ON MA.METHOD_AC = ME.METHOD_AC
          INNER JOIN INTERPRO.PROTEIN P ON MA.PROTEIN_AC = P.PROTEIN_AC
        UNION ALL
        SELECT
          ME.PROTEIN_AC,
          ME.CODE,
          ME.POS_FROM,
          ME.POS_TO,
          'T',
          'm',
          NULL,
          P.LEN,
          P.FRAGMENT
        FROM INTERPRO.MEROPS ME
          INNER JOIN INTERPRO.PROTEIN P ON ME.PROTEIN_AC = P.PROTEIN_AC
        ORDER BY PROTEIN_AC
        """
    )

    cur2 = con.cursor()  # second cursor for INSERT statements
    matches_all = []
    matches_filtered = []
    agg = {}
    cnt_proteins = 0
    cnt_matches = 0
    protein = None
    ts = time.time()
    for row in cur:
        protein_ac = row[0]

        if protein_ac != protein:
            # New protein! Flush previous protein
            if protein is not None:
                for method_ac in agg:
                    # Merge matches
                    min_pos = None
                    max_pos = None

                    for pos_start, pos_end in agg[method_ac]:
                        if min_pos is None or pos_start < min_pos:
                            min_pos = pos_start
                        if max_pos is None or pos_end > max_pos:
                            max_pos = pos_end

                    matches_filtered.append((method_ac, min_pos, max_pos))

                if matches_filtered:
                    proteins.put((protein, matches_filtered))

            matches_filtered = []
            agg = {}
            protein = protein_ac
            cnt_proteins += 1
            if not cnt_proteins % 1000000:
                logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt_proteins, cnt_proteins // (time.time() - ts)))

        method_ac = row[1]
        pos_start = math.floor(row[2])
        pos_end = math.floor(row[3])
        status = row[4]
        dbcode = row[5]
        method_type = row[6]
        protein_length = math.floor(row[7])
        is_fragment = row[8] == 'Y'

        matches_all.append((protein, method_ac, pos_start, pos_end, status, dbcode))

        if len(matches_all) == max_size:
            cnt_matches += max_size
            for i in range(0, len(matches_all), chunk_size):
                cur2.executemany(
                    """
                    INSERT /*+APPEND*/ INTO {}.MATCH (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, STATUS, DBCODE)
                    VALUES (:1, :2, :3, :4, :5, :6)
                    """.format(schema),
                    matches_all[i:i + chunk_size]
                )
                con.commit()
            matches_all = []

        if status != 'T' or is_fragment or dbcode == 'g':
            continue
        elif dbcode not in ('F', 'V'):
            matches_filtered.append((method_ac, pos_start, pos_end))
        elif method_ac not in agg:
            if method_type == 'F':
                # Families: use the entire protein sequence
                agg[method_ac] = [(1, protein_length)]
            else:
                agg[method_ac] = [(pos_start, pos_end)]
        elif method_type != 'F':
            # Since families use the entire protein, we can skip their matches
            agg[method_ac].append((pos_start, pos_end))

    cur.close()  # close the SELECT cursor as we don't need it anymore

    # Flush last protein
    for method_ac in agg:
        # Merge matches
        min_pos = None
        max_pos = None

        for pos_start, pos_end in agg[method_ac]:
            if min_pos is None or pos_start < min_pos:
                min_pos = pos_start
            if max_pos is None or pos_end > max_pos:
                max_pos = pos_end

        matches_filtered.append((protein, method_ac, min_pos, max_pos))

    if matches_filtered:
        proteins.put((protein, matches_filtered))
    matches_filtered = []
    agg = {}

    for i in range(0, len(matches_all), chunk_size):
        cur2.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.MATCH (PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, STATUS, DBCODE)
            VALUES (:1, :2, :3, :4, :5, :6)
            """.format(schema),
            matches_all[i:i + chunk_size]
        )
        con.commit()
    cnt_matches += len(matches_all)
    matches_all = []

    cnt_proteins += 1
    logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt_proteins, cnt_proteins // (time.time() - ts)))

    # Wait for comparators to finish
    for _ in comparators:
        proteins.put(None)

    for p in comparators:
        p.join()

    logging.info('{} matches inserted'.format(cnt_matches))
    comparisons.put(None)

    logging.info('indexing MATCH table')
    cur2.execute('CREATE INDEX I_MATCH$PROTEIN ON {}.MATCH (PROTEIN_AC) NOLOGGING'.format(schema))
    cur2.execute('CREATE INDEX I_MATCH$DBCODE ON {}.MATCH (DBCODE) NOLOGGING'.format(schema))

    logging.info('optimizing MATCH table')
    oracledb.gather_stats(cur2, schema, 'MATCH', cascade=True)
    oracledb.grant(cur2, 'SELECT', schema, 'MATCH', 'INTERPRO_SELECT')
    cur2.close()
    con.close()
    logging.info('MATCH table ready')

    aggregator.join()
    logging.info('complete')


def report_swiss_descriptions(dsn, schema, **kwargs):
    filepath = kwargs.get('output')

    if not isinstance(filepath, str) or not filepath:
        logging.critical("invalid 'output' argument")
        exit(1)

    dirname = os.path.dirname(os.path.abspath(filepath))
    if not os.path.dirname(dirname):
        os.makedirs(dirname)

    con = oracledb.connect(dsn)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, D.TEXT
          FROM {0}.METHOD2PROTEIN MP
        INNER JOIN {0}.METHOD M ON MP.METHOD_AC = M.METHOD_AC
        INNER JOIN {0}.DESC_VALUE D ON MP.DESC_ID = D.DESC_ID
        INNER JOIN {0}.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        WHERE MP.DBCODE = 'S'  
        """.format(schema)
    )

    old_entries = {}
    for acc, descr in cur:
        if acc in old_entries:
            old_entries[acc].append(descr)
        else:
            old_entries[acc] = [descr]

    cur.close()
    con.close()

    con = oracledb.connect(dsn)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, M.DESCRIPTION
        FROM {0}.METHOD2SWISS_DE M 
        INNER JOIN {0}.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        """.format(schema)
    )

    new_entries = {}
    for acc, descr in cur:
        if acc in new_entries:
            new_entries[acc].append(descr)
        else:
            new_entries[acc] = [descr]

    cur.execute('SELECT ENTRY_AC, ENTRY_TYPE, CHECKED FROM INTERPRO_ANALYSIS.ENTRY')
    entry_info = {acc: (entry_type, is_checked) for acc, entry_type, is_checked in cur}

    cur.close()
    con.close()

    entries = {}
    for acc in old_entries:
        try:
            entry_type, is_checked = entry_info[acc]
        except KeyError:
            continue

        if acc in new_entries:
            new_entry_descrs = new_entries.pop(acc)
            entries[acc] = {
                'acc': acc,
                'type': entry_type,
                'checked': is_checked,
                'lost': [descr for descr in old_entries[acc] if descr not in new_entry_descrs],
                'gained': [descr for descr in new_entries[acc] if descr not in new_entry_descrs]
            }
        else:
            # All descriptions lost (as we lost the entry)
            entries[acc] = {
                'acc': acc,
                'type': entry_type,
                'checked': is_checked,
                'lost': old_entries[acc],
                'gained': []
            }

    for acc in new_entries:
        entry_type, is_checked = entry_info[acc]
        entries[acc] = {
            'acc': acc,
            'type': entry_type,
            'checked': is_checked,
            'lost': [],
            'gained': new_entries[acc]
        }

    with open(filepath, 'wt') as fh:
        fh.write('Entry\tType\tChecked\t#Gained\t#Lost\tGained\tLost\n')

        for entry in sorted(entries.values(), key=lambda e: (0 if e['type'] == 'F' else 1, e['type'], e['acc'])):
            fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                entry['acc'],
                entry['type'],
                entry['checked'],
                len(entry['gained']),
                len(entry['lost']),
                ' | '.join(entry['gained']),
                ' | '.join(entry['lost'])
            ))
