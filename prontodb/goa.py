#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .oracledb import Connection


def load_annotations(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "PROTEIN2GO")

    con.execute(
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
    con.execute(
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

    con.execute(
        """
        CREATE INDEX I_PROTEIN2GO$P$G 
        ON {}.PROTEIN2GO (PROTEIN_AC, GO_ID) NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_PROTEIN2GO$E 
        ON {}.PROTEIN2GO (EVIDENCE) NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_PROTEIN2GO$RC 
        ON {}.PROTEIN2GO (REF_DB_CODE) NOLOGGING
        """.format(schema)
    )

    con.optimize_table(schema, "PROTEIN2GO", cascade=True)
    con.grant("SELECT", schema, "PROTEIN2GO", "INTERPRO_SELECT")


def load_publications(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "PUBLICATION")
    con.execute(
        """
        CREATE TABLE {}.PUBLICATION
        (
            ID VARCHAR2(25) NOT NULL,
            TITLE VARCHAR2(1500),
            FIRST_PUBLISHED_DATE DATE
        ) NOLOGGING        
        """.format(schema)
    )

    con.execute(
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

    con.execute(
        """
        ALTER TABLE {}.PUBLICATION
        ADD CONSTRAINT PK_PUBLICATION PRIMARY KEY (ID)
        """.format(schema)
    )

    con.optimize_table(schema, "PUBLICATION", cascade=True)
    con.grant("SELECT", schema, "PUBLICATION", "INTERPRO_SELECT")


def load_terms(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "TERM")
    con.execute(
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

    con.execute(
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

    con.execute(
        """
        ALTER TABLE {}.TERM
        ADD CONSTRAINT PK_TERM PRIMARY KEY (GO_ID)        
        """.format(schema)
    )

    con.optimize_table(schema, "TERM", cascade=True)
    con.grant("SELECT", schema, "TERM", "INTERPRO_SELECT")