#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import math
import os
from multiprocessing import Process, Queue
from tempfile import mkstemp

from .oracledb import Connection


def create_synonyms(dsn, src, dst, tables=[]):
    con = Connection(dsn)
    for table_name in tables:
        query = ("CREATE OR REPLACE SYNONYM {0}.{2} "
                 "FOR {1}.{2}".format(dst, src, table_name))
        con.execute(query)


def load_databases(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "CV_DATABASE")
    con.execute(
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

    con.execute(
        """
        INSERT /*+APPEND*/ INTO {}.CV_DATABASE (DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE)
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        """.format(schema)
    )
    con.commit()

    con.optimize_table(schema, "CV_DATABASE", cascade=True)
    con.grant("SELECT", schema, "CV_DATABASE", "INTERPRO_SELECT")


def load_signatures(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "METHODS")
    con.execute(
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

    con.execute(
        """
        INSERT /*+APPEND*/ INTO {}.METHOD (METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE)
        SELECT METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE
        FROM INTERPRO.METHOD
        """.format(schema)
    )
    con.commit()

    con.execute(
        """
        ALTER TABLE {}.METHOD
        ADD CONSTRAINT PK_METHOD PRIMARY KEY (METHOD_AC)
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_METHOD$DBCODE 
        ON {}.METHOD (DBCODE) NOLOGGING
        """.format(schema)
    )

    con.optimize_table(schema, "METHOD", cascade=True)
    con.grant("SELECT", schema, "METHOD", "INTERPRO_SELECT")


def load_taxa(dsn, schema, chunk_size=100000):
    con = Connection(dsn)
    con.drop_table(schema, "ETAXI")
    con.drop_table(schema, "LINEAGE")
    con.execute(
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

    con.execute(
        """
        CREATE TABLE {}.LINEAGE
        (
            LEFT_NUMBER NUMBER NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            RANK VARCHAR2(50)
        ) NOLOGGING
        """.format(schema)
    )

    con.execute(
        """
        INSERT /*+APPEND*/ INTO {}.ETAXI (
            TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, 
            LEFT_NUMBER, RIGHT_NUMBER, FULL_NAME
        )
        SELECT 
            TAX_ID, PARENT_ID, SCIENTIFIC_NAME, RANK, 
            LEFT_NUMBER, RIGHT_NUMBER, FULL_NAME
        FROM INTERPRO.ETAXI
        """.format(schema)
    )
    con.commit()
    con.execute(
        """
        ALTER TABLE {}.ETAXI
        ADD CONSTRAINT PK_ETAXI PRIMARY KEY (TAX_ID)
        """.format(schema)
    )

    con.optimize_table(schema, "ETAXI", cascade=True)
    con.grant("SELECT", schema, "ETAXI", "INTERPRO_SELECT")

    taxons = {}
    query = """
        SELECT TAX_ID, PARENT_ID, LEFT_NUMBER, RANK
        FROM {}.ETAXI
    """.format(schema)
    for tax_id, parent_id, left_num, rank in con.get(query):
        if tax_id != 131567:
            """
            taxID 131567 (cellular organisms) contains three superkingdoms:
                * Bacteria (2)
                * Archaea (2157)
                * Eukaryota (2759)

            therefore it is not needed (we don't want a meta-superkingdom)
            """
            taxons[tax_id] = {
                "parent": parent_id,
                "left_num": left_num,
                "rank": "superkingdom" if parent_id == 1 else rank
            }

    lineage = []
    for tax_id in taxons:
        t = taxons[tax_id]
        left_num = t["left_num"]
        if not left_num:
            continue
        elif t["rank"] != "no rank":
            lineage.append((left_num, tax_id, t["rank"]))

        parent_id = t["parent"]
        while parent_id:
            if parent_id not in taxons:
                # taxID 131567 missing from dictionary
                break

            t = taxons[parent_id]
            if t["rank"] != "no rank":
                lineage.append((left_num, parent_id, t["rank"]))

            parent_id = t["parent"]

    for i in range(0, len(lineage), chunk_size):
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.LINEAGE (LEFT_NUMBER, TAX_ID, RANK)
            VALUES (:1, :2, :3)
            """.format(schema),
            lineage[i:i+chunk_size]
        )
        con.commit()

    con.execute(
        """
        CREATE INDEX I_LINEAGE$L$R
        ON {}.LINEAGE (LEFT_NUMBER, RANK)
        NOLOGGING
        """.format(schema)
    )

    con.optimize_table(schema, "LINEAGE", cascade=True)
    con.grant("SELECT", schema, "LINEAGE", "INTERPRO_SELECT")


def load_proteins(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "PROTEIN")
    con.execute(
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

    con.execute(
        """
        INSERT /*+APPEND*/ INTO {}.PROTEIN (
            PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, TAX_ID
        )
        SELECT 
            PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, TAX_ID
        FROM INTERPRO.PROTEIN
        """.format(schema)
    )
    con.commit()

    con.execute(
        """
        ALTER TABLE {}.PROTEIN
        ADD CONSTRAINT PK_PROTEIN PRIMARY KEY (PROTEIN_AC)
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_PROTEIN$DBCODE 
        ON {}.PROTEIN (DBCODE) NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_PROTEIN$NAME 
        ON {}.PROTEIN (NAME) NOLOGGING
        """.format(schema)
    )

    con.optimize_table(schema, "PROTEIN", cascade=True)
    con.grant("SELECT", schema, "PROTEIN", "INTERPRO_SELECT")


class MatchComparator(Process):
    def __init__(self, tasks, results, max_gap):
        super().__init__()
        self.tasks = tasks
        self.results = results
        self.max_gap = max_gap

    def run(self):
        while True:
            task = self.tasks.get()
            if task is None:
                break

            acc, dbcode, length, descr_id, left_num, matches = task
            structure = self.condense(matches, self.max_gap)
            signatures, comparisons = self.compare(matches)
            self.results.put((
                acc, dbcode, length, descr_id, left_num,
                structure, signatures, comparisons
            ))

    @staticmethod
    def compare(matches):
        comparisons = []
        signatures = {}
        for m1_acc, m1_start, m1_end in matches:
            m1_len = m1_end - m1_start + 1

            if m1_acc in signatures:
                signatures[m1_acc] += 1
            else:
                signatures[m1_acc] = 1

            for m2_acc, m2_start, m2_end in matches:
                m2_len = m2_end - m2_start + 1
                """
                1           13      end position *is* included
                -------------       match 1 (13 - 1 + 1 = 13 aa)
                        --------    match 2 (16 - 9 + 1 = 8 aa)
                        9      16
                                    overlap = 13 - 9 + 1 = 5
                                    frac 1 = 5 / 13 = 0.38...
                                    frac 2 = 5 / 8 = 0.625
                """
                if m1_acc <= m2_acc:
                    comparisons.append((
                        m1_acc,
                        m1_len,
                        m2_acc,
                        m2_len,
                        # overlap length
                        min(m1_end, m2_end) - max(m1_start, m2_start) + 1
                    ))

        return list(signatures.items()), comparisons

    @staticmethod
    def condense(matches, max_gap):
        # flatten matches
        locations = []

        for method_ac, pos_start, pos_end in matches:
            locations.append((pos_start, method_ac))
            locations.append((pos_end, method_ac))

        """
        Evaluate the protein's match structure, 
            i.e. how signatures match the proteins

        -----------------------------   Protein
         <    >                         Signature 1
           <    >                       Signature 2
                      < >               Signature 3

        Flattened:
        -----------------------------   Protein
         < <  > >     < >
         1 2  1 2     3 3

        Structure, with '-' representing a "gap" 
            (more than N bp between two positions):
        1212-33
        """

        # Sort locations by position
        locations.sort()

        """
        Do not set the offset to 0, but to the first position:
        if two proteins have the same structure,
        but the first position of one protein is > max_gap 
        while the first position of the other protein is <= max_gap,
        a gap will be used for the first protein and not for the other, 
        which will results in two different structures
        """
        offset = locations[0][0]

        structure = []  # overall match structure
        methods = []  # close signatures (less than max_gap between two positions)

        for pos, method_ac in sorted(locations):
            if pos > offset + max_gap:
                for _pos, _ac in methods:
                    structure.append(_ac)

                structure.append('')  # add a gap
                methods = []

            offset = pos
            methods.append((pos, method_ac))

        for _pos, _ac in methods:
            structure.append(_ac)

        return '/'.join(structure)


class MatchAggregator(Process):
    def __init__(self, dsn, schema, tasks, tmpdir=None, chunk_size=100000):
        super().__init__()
        self.dsn = dsn
        self.schema = schema
        self.tasks = tasks
        self.tmpdir = tmpdir
        self.chunk_size = chunk_size

    def run(self):
        structures = {}
        signatures = {}
        comparisons = {}
        proteins = []
        files = []

        while True:
            task = self.tasks.get()
            if task is None:
                break

            (accession, dbcode, length, descr_id, left_num,
             _structure, _signatures, _comparisons) = task

            # struct: condensed match structure of the protein
            if _structure in structures:
                code = structures[_structure]
            else:
                code = self.base36encode(len(structures) + 1)
                structures[_structure] = code

            p = {
                "acc": accession,
                "dbcode": dbcode,
                "length": length,
                "descr": descr_id,
                "leftnum": left_num,
                "code": code,
                "signatures": []
            }
            proteins.append(p)

            # _signatures: number of times each signature matched the protein
            for acc, n_matches in _signatures:
                p["signatures"].append(acc)

                if acc not in signatures:
                    signatures[acc] = {
                        "proteins": 0,
                        "matches": 0
                    }

                signatures[acc]['proteins'] += 1
                signatures[acc]['matches'] += n_matches

            if len(proteins) == self.chunk_size:
                files.append(self.dump(proteins, self.tmpdir))
                proteins = []

            # _comparisons: match overlaps between signatures
            collocations = set()
            overlaps = set()

            for acc_1, len_1, acc_2, len_2, overlap in _comparisons:
                if acc_1 in comparisons:
                    d = comparisons[acc_1]
                else:
                    d = comparisons[acc_1] = {}

                if acc_2 in d:
                    comp = d[acc_2]
                else:
                    comp = d[acc_2] = {
                        # number of proteins in which both signatures occur (not necessarily overlapping)
                        'prot': 0,
                        # number of proteins in which signatures overlap
                        'prot_over': 0,
                        # number of times signatures overlap (>= prot_over)
                        'over': 0,
                        # sum of overlap lengths (to compute the average length later)
                        'length': 0,
                        # sum of fractions of matches overlapping (overlap length / match length)
                        'frac_1': 0,
                        'frac_2': 0
                    }

                acc = (acc_1, acc_2)
                if acc not in collocations:
                    # First collocation (acc_1, acc_2) for the current protein
                    collocations.add(acc)
                    comp['prot'] += 1

                if overlap > (min(len_1, len_2) / 2):
                    """
                    Consider that matches significantly overlap 
                    if the overlap is longer 
                    than the half of the shortest match
                    """
                    comp['frac_1'] += overlap / len_1
                    comp['frac_2'] += overlap / len_2
                    comp['over'] += 1
                    comp['length'] += overlap

                    if acc not in overlaps:
                        # First overlap (acc_1, acc_2) for the current protein
                        overlaps.add(acc)
                        comp['prot_over'] += 1

        # Dump proteins left in memory
        if proteins:
            files.append(self.dump(proteins))
            proteins = []

        con = Connection(self.dsn)
        candidates = set()
        non_prosite_candidates = set()
        query = """
            SELECT METHOD_AC, DBCODE
            FROM INTERPRO.METHOD
            WHERE CANDIDATE = 'Y'
        """
        for accession, dbcode in con.get(query):
            candidates.add(accession)
            if dbcode != 'P':
                non_prosite_candidates.add(accession)

        # Determine relations/adjacent
        relations = {}
        adjacents = {}
        for acc_1 in comparisons:
            s_1 = signatures[acc_1]

            for acc_2 in comparisons[acc_1]:
                s_2 = signatures[acc_2]
                c = comparisons[acc_1][acc_2]

                if c['prot_over'] >= s_1['proteins'] * 0.4 or c['prot_over'] >= s_2['proteins'] * 0.4:
                    # is a relation
                    if acc_1 in relations:
                        relations[acc_1].append(acc_2)
                    else:
                        relations[acc_1] = [acc_2]
                elif c['prot'] >= s_1['proteins'] * 0.1:
                    # is adjacent
                    if acc_1 in adjacents:
                        adjacents[acc_1].append(acc_2)
                    else:
                        adjacents[acc_1] = [acc_2]

                if acc_1 != acc_2:
                    if c['prot_over'] >= s_1['proteins'] * 0.4 or c['prot_over'] >= s_2['proteins'] * 0.4:
                        # is a relation
                        if acc_2 in relations:
                            relations[acc_2].append(acc_1)
                        else:
                            relations[acc_2] = [acc_1]
                    elif c['prot'] >= s_2['proteins'] * 0.1:
                        # is adjacent
                        if acc_2 in adjacents:
                            adjacents[acc_2].append(acc_1)
                        else:
                            adjacents[acc_2] = [acc_1]

        """
        Determine extra/adjacent relation:

        Let A, and B be two signatures in a relationship.
        extra relation:     number of signatures in a relationship with A but not with B
        adjacent relation:  number of signatures in a relationship with A and adjacent to B
        """
        extra_relations = {}
        adj_relations = {}
        for acc_1 in relations:
            # signatures in a relationship with acc_1 (that are candidate and not from PROSITE pattern)
            rel_1 = set(relations[acc_1]) & non_prosite_candidates

            for acc_2 in relations[acc_1]:
                # signatures in a relationship with acc_2
                rel_2 = set(relations.get(acc_2, []))

                # signature adjacent to acc_2
                adj_2 = set(adjacents.get(acc_2, []))

                extra = rel_1 - rel_2
                if acc_1 in extra_relations:
                    extra_relations[acc_1][acc_2] = len(extra)
                else:
                    extra_relations[acc_1] = {acc_2: len(extra)}

                adj = rel_1 & adj_2
                if acc_1 in adj_relations:
                    adj_relations[acc_1][acc_2] = len(adj)
                else:
                    adj_relations[acc_1] = {acc_2: len(adj)}

        # Make predictions
        predictions = []
        for acc_1 in comparisons:
            s_1 = signatures[acc_1]

            if acc_1 not in candidates:
                continue

            for acc_2 in comparisons[acc_1]:
                if acc_1 == acc_2 or acc_2 not in candidates:
                    continue

                c = comparisons[acc_1][acc_2]
                s_2 = signatures[acc_2]

                if c['prot_over'] >= s_1['proteins'] * 0.4 or c['prot_over'] >= s_2['proteins'] * 0.4:
                    extra_1 = extra_relations.get(acc_1, {}).get(acc_2, 0)
                    extra_2 = extra_relations.get(acc_2, {}).get(acc_1, 0)

                    # inverting acc2 and acc1 to have the same predictions than HH, todo: check if it is really ok
                    adj_1 = adj_relations.get(acc_2, {}).get(acc_1, 0)
                    adj_2 = adj_relations.get(acc_1, {}).get(acc_2, 0)

                    if c['over']:
                        """
                        frac_1 and frac_2 are the sum of ratios (overlap / match length).
                        if close to 1, it means that the match was mostly contained by the overlap in most cases
                                A   -------------
                                B       ---             

                            frac(B) = overlap(A,B) / length(B) = 1 (indeed B is 100% within the overlap)

                        len_1 and len_2 are the average of sums.
                        If len(B) is close to 1, it means that B was mostly within the overlap in most cases
                            so B < A (because B ~ overlap and A >= overlap)
                            so B CONTAINED_BY A 
                        """
                        len_1 = c['frac_1'] / c['over']
                        len_2 = c['frac_2'] / c['over']
                    else:
                        len_1 = len_2 = 0

                    """
                    Parent/Child relationships:
                    The protein/matches made by the child entry must be a complete (>75%) subset of the parent entry
                    if over(A) > 0.75, it means that A overlaps with B in at least 75% of its proteins or matches:
                        A is a CHILD_OF B
                    """
                    over_1 = min(c['over'] / s_1['matches'], c['prot_over'] / s_1['proteins'])
                    over_2 = min(c['over'] / s_2['matches'], c['prot_over'] / s_2['proteins'])

                    if len_1 >= 0.5 and len_2 >= 0.5:
                        if over_1 > 0.75 and over_2 >= 0.75 and not extra_1 and not extra_2:
                            prediction = 'ADD_TO'
                        elif over_1 > 0.75 and not extra_1 and not adj_1:
                            prediction = 'CHILD_OF'
                        elif over_2 > 0.75 and not extra_2 and not adj_2:
                            prediction = 'PARENT_OF'  # acc2 child of acc1
                        elif len_1 >= 0.9:
                            prediction = 'C/C' if len_2 >= 0.9 else 'CONTAINED_BY'
                        elif len_2 >= 0.9:
                            prediction = 'CONTAINER_OF'
                        else:
                            prediction = 'OVERLAPS'
                    elif len_1 >= 0.9:
                        prediction = 'CONTAINED_BY'
                    elif len_2 >= 0.9:
                        prediction = 'CONTAINER_OF'
                    else:
                        prediction = 'OVERLAPS'

                    predictions.append((acc_1, acc_2, prediction))

                    # switch (acc_1, acc_2) -> (acc_2, acc_1)
                    if prediction == 'CHILD_OF':
                        prediction = 'PARENT_OF'
                    elif prediction == 'PARENT_OF':
                        prediction = 'CHILD_OF'
                    elif prediction == 'CONTAINED_BY':
                        prediction = 'CONTAINER_OF'
                    elif prediction == 'CONTAINER_OF':
                        prediction = 'CONTAINED_BY'
                    predictions.append((acc_2, acc_1, prediction))

        # Populating METHOD_MATCH
        con.drop_table(self.schema, "METHOD_MATCH")
        con.execute(
            """
            CREATE TABLE {}.METHOD_MATCH
            (
                METHOD_AC VARCHAR2(25) NOT NULL,
                N_MATCHES NUMBER NOT NULL,
                N_PROT NUMBER NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        signatures = [(acc, s['matches'], s['proteins']) for acc, s in signatures.items()]
        for i in range(0, len(signatures), self.chunk_size):
            con.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.METHOD_MATCH (METHOD_AC, N_MATCHES, N_PROT)
                VALUES (:1, :2, :3)
                """.format(self.schema),
                signatures[i:i + self.chunk_size]
            )
            con.commit()

        # Optimizing METHOD_MATCH
        con.execute(
            """
            ALTER TABLE {}.METHOD_MATCH
            ADD CONSTRAINT PK_METHOD_MATCH PRIMARY KEY (METHOD_AC)
            """.format(self.schema)
        )
        con.optimize_table(self.schema, "METHOD_MATCH", cascade=True)
        con.grant("SELECT", self.schema, "METHOD_MATCH", "INTERPRO_SELECT")

        # Creating METHOD_OVERLAP
        con.drop_table(self.schema, "METHOD_OVERLAP")
        con.execute(
            """
            CREATE TABLE {}.METHOD_OVERLAP
            (
                METHOD_AC1 VARCHAR2(25) NOT NULL,
                METHOD_AC2 VARCHAR2(25) NOT NULL,
                N_PROT NUMBER NOT NULL,
                N_OVER NUMBER NOT NULL,
                N_PROT_OVER NUMBER NOT NULL,
                AVG_OVER NUMBER NOT NULL,
                AVG_FRAC1 NUMBER NOT NULL,
                AVG_FRAC2 NUMBER NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        # Computing data for METHOD_OVERLAP
        overlaps = []
        for acc_1 in comparisons:
            for acc_2 in comparisons[acc_1]:
                c = comparisons[acc_1][acc_2]

                if c['over']:
                    avg_frac1 = 100 * c['frac_1'] / c['over']
                    avg_frac2 = 100 * c['frac_2'] / c['over']
                    avg_over = c['length'] / c['over']
                else:
                    avg_frac1 = avg_frac2 = avg_over = 0

                """
                Cast to float as cx_Oracle 6.1 throws TypeError (expecting integer)
                when the value of the 1st record is an integer (e.g. AVG_OVER=0)
                and the value for the second is not (e.g. AVG_OVER=25.6)
                """
                overlaps.append((
                    acc_1, acc_2, c['prot'], c['over'], c['prot_over'],
                    float(avg_over), float(avg_frac1), float(avg_frac2)
                ))

                if acc_1 != acc_2:
                    overlaps.append((
                        acc_2, acc_1, c['prot'], c['over'], c['prot_over'],
                        float(avg_over), float(avg_frac2), float(avg_frac1)
                    ))

        # Populating METHOD_OVERLAP
        for i in range(0, len(overlaps), self.chunk_size):
            con.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.METHOD_OVERLAP (
                    METHOD_AC1, METHOD_AC2, N_PROT, N_OVER, 
                    N_PROT_OVER, AVG_OVER, AVG_FRAC1, AVG_FRAC2
                )
                VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
                """.format(self.schema),
                overlaps[i:i + self.chunk_size]
            )
            con.commit()

        # Optimizing METHOD_OVERLAP
        con.execute(
            """
            ALTER TABLE {}.METHOD_OVERLAP
            ADD CONSTRAINT PK_METHOD_OVERLAP 
            PRIMARY KEY (METHOD_AC1, METHOD_AC2)
            """.format(self.schema)
        )
        con.optimize_table(self.schema, "METHOD_OVERLAP", cascade=True)
        con.grant("SELECT", self.schema, "METHOD_OVERLAP", "INTERPRO_SELECT")

        # Creating METHOD_PREDICTION
        con.drop_table(self.schema, "METHOD_PREDICTION")
        con.execute(
            """
            CREATE TABLE {}.METHOD_PREDICTION
            (
                METHOD_AC1 VARCHAR2(25) NOT NULL,
                METHOD_AC2 VARCHAR2(25) NOT NULL,
                RELATION VARCHAR2(15) NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        # Populating METHOD_PREDICTION
        for i in range(0, len(predictions), self.chunk_size):
            con.executemany(
                """
                INSERT /*+APPEND*/ 
                INTO {}.METHOD_PREDICTION (METHOD_AC1, METHOD_AC2, RELATION)
                VALUES (:1, :2, :3)
                """.format(self.schema),
                predictions[i:i + self.chunk_size]
            )
            con.commit()

        # Optimizing METHOD_PREDICTION
        con.execute(
            """
            ALTER TABLE {}.METHOD_PREDICTION
            ADD CONSTRAINT PK_METHOD_PREDICTION 
            PRIMARY KEY (METHOD_AC1, METHOD_AC2)
            """.format(self.schema)
        )
        con.optimize_table(self.schema, "METHOD_PREDICTION", cascade=True)
        con.grant("SELECT", self.schema, "METHOD_PREDICTION", "INTERPRO_SELECT")

        # Get lineages
        ranks = {
            "superkingdom",
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species"
        }
        left_numbers = {}
        query = """
            SELECT RANK_LEFT_NUMBER_TAX_ID
            FROM {}.LINEAGE
        """.format(self.schema)
        for rank, left_num, tax_id in con.get(query):
            if rank not in ranks:
                continue
            elif left_num not in left_numbers:
                left_numbers[left_num] = {}
            left_numbers[left_num][rank] = tax_id

        # Creating METHOD2PROTEIN
        con.drop_table(self.schema, "METHOD2PROTEIN")
        con.execute(
            """
            CREATE TABLE {}.METHOD2PROTEIN
            (
                METHOD_AC VARCHAR2(25) NOT NULL,
                PROTEIN_AC VARCHAR2(15) NOT NULL,
                DBCODE CHAR(1) NOT NULL,
                CONDENSE VARCHAR(100) NOT NULL,
                LEN NUMBER(5) NOT NULL,
                LEFT_NUMBER NUMBER NOT NULL,
                DESC_ID NUMBER(10) NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        # Populating METHOD2PROTEIN_STG by loading files
        methods = {}
        for filepath in files:
            with gzip.open(filepath, 'rt') as fh:
                proteins = json.load(fh)

            os.unlink(filepath)

            data = []
            for p in proteins:
                left_num = p["leftnum"]
                items = list(left_numbers.get(left_num, {}).items())

                for method_ac in p["signatures"]:
                    data.append((
                        method_ac,
                        p["acc"],
                        p["dbcode"],
                        p["code"],
                        p["length"],
                        left_num,
                        p["descr"]
                    ))

                    # Counting the taxonomic origin of signatures
                    if method_ac in methods:
                        ranks = methods[method_ac]
                    else:
                        ranks = methods[method_ac] = {}

                    for rank, tax_id in items:
                        if rank not in ranks:
                            ranks[rank] = {tax_id: 1}
                        elif tax_id not in ranks[rank]:
                            ranks[rank][tax_id] = 1
                        else:
                            ranks[rank][tax_id] += 1

            for i in range(0, len(data), self.chunk_size):
                con.executemany(
                    """
                    INSERT /*+APPEND*/ INTO {}.METHOD2PROTEIN (
                        METHOD_AC, PROTEIN_AC, DBCODE, 
                        CONDENSE, LEN, LEFT_NUMBER, DESC_ID
                    )
                    VALUES (:1, :2, :3, :4, :5, :6, :7)
                    """.format(self.schema),
                    data[i:i + self.chunk_size]
                )
                con.commit()

            # Can I haz some memory back?
            proteins = None
            data = None

        # Optimizing METHOD2PROTEIN
        con.execute(
            """
            ALTER TABLE {}.METHOD2PROTEIN
            ADD CONSTRAINT PK_METHOD2PROTEIN 
            PRIMARY KEY (METHOD_AC, PROTEIN_AC)
            """.format(self.schema)
        )
        con.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$M
            ON {}.METHOD2PROTEIN (METHOD_AC)
            NOLOGGING
            """.format(self.schema)
        )
        con.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$P
            ON {}.METHOD2PROTEIN (PROTEIN_AC)
            NOLOGGING
            """.format(self.schema)
        )
        con.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$LN
            ON {}.METHOD2PROTEIN (LEFT_NUMBER)
            NOLOGGING
            """.format(self.schema)
        )
        con.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$M$DB
            ON {}.METHOD2PROTEIN (METHOD_AC, DBCODE)
            NOLOGGING
            """.format(self.schema)
        )

        con.optimize_table(self.schema, "METHOD2PROTEIN", cascade=True)
        con.grant("SELECT", self.schema, "METHOD2PROTEIN", "INTERPRO_SELECT")

        # Creating METHOD_TAXA
        con.drop_table(self.schema, "METHOD_TAXA")
        con.execute(
            """
            CREATE TABLE {}.METHOD_TAXA
            (
                METHOD_AC VARCHAR2(25) NOT NULL,
                RANK VARCHAR2(50) NOT NULL,
                TAX_ID NUMBER(10) NOT NULL,
                PROTEIN_COUNT NUMBER(10) NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        # Populating METHOD_TAXA
        for method_ac, ranks in methods.items():
            con.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.METHOD_TAXA (
                    METHOD_AC, RANK, TAX_ID, PROTEIN_COUNT
                )
                VALUES (:1, :2, :3, :4)
                """.format(self.schema),
                [
                    (method_ac, rank, tax_id, count)
                    for rank, taxa in ranks.items()
                    for tax_id, count in taxa.items()
                ]
            )
            con.commit()

        # Optimizing METHOD_TAXA
        con.execute(
            """
            CREATE INDEX I_METHOD_TAXA
            ON {}.METHOD_TAXA (METHOD_AC, RANK)
            NOLOGGING
            """.format(self.schema)
        )

        con.optimize_table(self.schema, "METHOD_TAXA", cascade=True)
        con.grant("SELECT", self.schema, "METHOD_TAXA", "INTERPRO_SELECT")

    @staticmethod
    def base36encode(num):
        chars = "0123456789abcdefghijklmnopqrstuvwxyz"
        encoded = ""

        while num > 0:
            num, remainder = divmod(num, 36)
            encoded = chars[remainder] + encoded

        return encoded

    @staticmethod
    def dump(data, tmpdir=None):
        fd, filepath = mkstemp(suffix='.json.gz', dir=tmpdir)
        os.close(fd)

        with gzip.open(filepath, 'wt') as fh:
            json.dump(data, fh)

        return filepath


def load_matches(dsn, schema, **kwargs):
    threads = kwargs.get("threads", 1)
    max_gap = kwargs.get("max_gap", 20)
    tmpdir = kwargs.get("tmpdir")
    chunk_size = kwargs.get("chunk_size", 100000)
    limit = kwargs.get("limit", 0)

    q_proteins = Queue(maxsize=chunk_size)
    q_matches = Queue(maxsize=chunk_size)

    # -2: main thread, aggregator
    comparators = [
        MatchComparator(q_proteins, q_matches, max_gap)
        for _ in range(max(1, threads-2))
    ]

    aggregator = MatchAggregator(dsn, schema, q_matches,
                                 tmpdir=tmpdir, chunk_size=chunk_size)

    for p in comparators:
        p.start()

    aggregator.start()

    con = Connection(dsn)
    con.drop_table(schema, "MATCH")
    con.execute(
        """
        CREATE TABLE {}.MATCH
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            MODEL_AC VARCHAR2(25) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            FRAGMENTS VARCHAR2(200) DEFAULT NULL
        ) NOLOGGING
        """.format(schema)
    )
    
    query = """
        SELECT M.*, P.LEN, P.FRAGMENT, P.DBCODE, PD.DESC_ID, NVL(E.LEFT_NUMBER, 0)
        FROM (
            SELECT
              MA.PROTEIN_AC AS PROTEIN_AC,
              MA.METHOD_AC AS METHOD_AC,
              MA.MODEL_AC AS MODEL_AC,
              MA.POS_FROM AS POS_FROM,
              MA.POS_TO AS POS_TO,
              MA.FRAGMENTS AS FRAGMENTS,
              MA.DBCODE AS M_DBCODE,
              ME.SIG_TYPE AS SIG_TYPE
            FROM INTERPRO.MATCH MA
              INNER JOIN INTERPRO.METHOD ME ON MA.METHOD_AC = ME.METHOD_AC
            UNION ALL
            SELECT
              FM.PROTEIN_AC AS PROTEIN_AC,
              FM.METHOD_AC AS METHOD_AC,
              FM.METHOD_AC AS MODEL_AC,
              FM.POS_FROM AS POS_FROM,
              FM.POS_TO AS POS_TO,
              NULL AS FRAGMENTS,
              FM.DBCODE AS M_DBCODE,
              ME.SIG_TYPE AS SIG_TYPE
            FROM INTERPRO.FEATURE_MATCH FM
              INNER JOIN INTERPRO.METHOD ME ON FM.METHOD_AC = ME.METHOD_AC
            WHERE FM.DBCODE = 'g'
            UNION ALL            
            SELECT
              ME.PROTEIN_AC AS PROTEIN_AC,
              ME.CODE AS METHOD_AC,
              ME.CODE AS MODEL_AC,
              ME.POS_FROM AS POS_FROM,
              ME.POS_TO AS POS_TO,
              NULL AS FRAGMENTS,
              'm' AS M_DBCODE,
              NULL AS SIG_TYPE
            FROM INTERPRO.MEROPS ME        
        ) M
        INNER JOIN INTERPRO.PROTEIN P ON M.PROTEIN_AC = P.PROTEIN_AC
        INNER JOIN INTERPRO.ETAXI E ON P.TAX_ID = E.TAX_ID
        INNER JOIN {}.PROTEIN_DESC PD ON M.PROTEIN_AC = PD.PROTEIN_AC
        ORDER BY M.PROTEIN_AC
        """.format(schema)

    matches = []            # matches to be inserted in the MATCH table
    matches_predict = []    # matches to be used for protein match structure and predictions
    methods = {}            # methods having matches to merge
    protein = None          # previous protein accession
    prot_dbcode = None      # protein DB code (S: SwissProt, T: TrEMBL)
    length = None           # Sequence length
    descr_id = None         # Description ID
    left_num = None         # Taxon left number
    cnt = 0
    for row in con.get(query):
        protein_acc = row[0]

        if protein_acc != protein:
            if protein is not None:
                for method_acc in methods:
                    # Merge matches
                    min_pos = None
                    max_pos = None

                    for start, end in methods[method_acc]:
                        if min_pos is None or start < min_pos:
                            min_pos = start
                        if max_pos is None or end > max_pos:
                            max_pos = end

                    matches_predict.append((method_acc, min_pos, max_pos))

                if matches_predict:
                    q_proteins.put((protein, prot_dbcode, length, descr_id,
                                    left_num, matches_predict))

                cnt += 1
                if cnt == limit:
                    break

            matches_predict = []
            methods = {}
            protein = protein_acc

        method_acc = row[1]
        model_acc = row[1] if row[2] is None else row[2]
        start = math.floor(row[3])
        end = math.floor(row[4])
        fragments = row[5]
        method_dbcode = row[6]
        method_type = row[7]
        length = math.floor(row[8])
        is_fragment = row[9] == 'Y'
        prot_dbcode = row[10]
        descr_id = row[11]
        left_num = row[12]

        matches.append((
            protein, method_acc, model_acc, start, end,
            method_dbcode, fragments
        ))
        if len(matches) == chunk_size:
            con.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.MATCH (
                    PROTEIN_AC, METHOD_AC, MODEL_AC, 
                    POS_FROM, POS_TO, DBCODE, FRAGMENTS
                )
                VALUES (:1, :2, :3, :4, :5, :6, :7)
                """.format(schema),
                matches
            )
            con.commit()
            matches = []

        """
        MobiDB-lite: keep matches for MATCH table, but not for predictions

        PANTHER & PRINTS:
            Merge protein matches.
            If the signature is a family*, use the entire protein.

            * all PANTHER signatures, almost all PRINTS signatures
        """
        if is_fragment or method_dbcode == 'g':
            continue
        elif method_dbcode not in ('F', 'V'):
            matches_predict.append((method_acc, start, end))
        elif method_acc not in methods:
            if method_type == 'F':
                # Families: use the entire protein sequence
                methods[method_acc] = [(1, length)]
            else:
                methods[method_acc] = [(start, end)]
        elif method_type != 'F':
            # Since families use the entire protein, skip their matches
            methods[method_acc].append((start, end))

    # Flush last protein
    for method_acc in methods:
        # Merge matches
        min_pos = None
        max_pos = None

        for start, end in methods[method_acc]:
            if min_pos is None or start < min_pos:
                min_pos = start
            if max_pos is None or end > max_pos:
                max_pos = end

        matches_predict.append((method_acc, min_pos, max_pos))

    if matches_predict:
        q_proteins.put((protein, prot_dbcode, length, descr_id,
                        left_num, matches_predict))

    if matches:
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.MATCH (
                PROTEIN_AC, METHOD_AC, MODEL_AC, 
                POS_FROM, POS_TO, DBCODE, FRAGMENTS
            )
            VALUES (:1, :2, :3, :4, :5, :6, :7)
            """.format(schema),
            matches
        )
        con.commit()
        matches = []

    for _ in comparators:
        q_proteins.put(None)
        
    # Wait for comparators to complete
    for p in comparators:
        p.join()

    # Triggers prediction/ METHOD2PROTEIN creating
    q_matches.put(None)

    # In the meantime: index MATCH
    con.execute(
        """
        CREATE INDEX I_MATCH$PROTEIN 
        ON {}.MATCH (PROTEIN_AC) NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_MATCH$DBCODE
        ON {}.MATCH (DBCODE) NOLOGGING
        """.format(schema)
    )
    con.optimize_table(schema, "MATCH", cascade=True)
    con.grant("SELECT", schema, "MATCH", "INTERPRO_SELECT")

    aggregator.join()


def copy_schema(dsn, schema):
    proc = "{}.copy_interpro_analysis.refresh".format(schema)
    Connection(dsn).exec(proc)


def clear_schema(dsn, schema):
    con = Connection(dsn)
    for table in con.get_tables(schema):
        con.drop_table(schema, table)


def report_description_changes(dsn, schema, output):
    con = Connection(dsn)

    """
    Get descriptions for before update 
    (METHOD2SWISS_DE was populated during protein update)
    """
    query = """
        SELECT DISTINCT EM.ENTRY_AC, M.DESCRIPTION
        FROM {0}.METHOD2SWISS_DE M
        INNER JOIN {0}.ENTRY2METHOD EM 
        ON M.METHOD_AC = EM.METHOD_AC
    """.format(schema)
    old_entries = {}
    for acc, text in con.get(query):
        if acc in old_entries:
            old_entries[acc].append(text)
        else:
            old_entries[acc] = [text]

    # Get descriptions for after update
    query = """
        SELECT DISTINCT EM.ENTRY_AC, D.TEXT
        FROM {0}.METHOD2PROTEIN MP
        INNER JOIN {0}.METHOD M ON MP.METHOD_AC = M.METHOD_AC
        INNER JOIN {0}.DESC_VALUE D ON MP.DESC_ID = D.DESC_ID
        INNER JOIN {0}.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        WHERE MP.DBCODE = 'S' 
    """.format(schema)
    new_entries = {}
    for acc, text in con.get(query):
        if acc in new_entries:
            new_entries[acc].append(text)
        else:
            new_entries[acc] = [text]

    query = """
        SELECT ENTRY_AC, ENTRY_TYPE, CHECKED
        FROM {}.ENTRY
    """.format(schema)
    entry_info = {
        acc: (entry_type, is_checked)
        for acc, entry_type, is_checked in con.get(query)
    }

    entries = {}
    for acc in old_entries:
        try:
            entry_type, is_checked = entry_info[acc]
        except KeyError:
            continue

        if acc in new_entries:
            descriptions = new_entries.pop(acc)
            entries[acc] = {
                "acc": acc,
                "type": entry_type,
                "checked": is_checked,
                "lost": [d for d in old_entries[acc]
                         if d not in descriptions],
                "gained": [d for d in descriptions
                           if d not in old_entries[acc]]
            }
        else:
            # All descriptions lost (we lost the entry)
            entries[acc] = {
                "acc": acc,
                "type": entry_type,
                "checked": is_checked,
                "lost": old_entries[acc],
                "gained": []
            }

    for acc in new_entries:
        entry_type, is_checked = entry_info[acc]
        entries[acc] = {
            "acc": acc,
            "type": entry_type,
            "checked": is_checked,
            "lost": [],
            "gained": new_entries[acc]
        }

    with open(output, "wt") as fh:
        fh.write("Entry\tType\tChecked\t#Gained\t#Lost\tGained\tLost\n")

        for e in sorted(entries.values(), key=lambda x: (0 if x["type"] == 'F' else 1, x["type"], x["acc"])):
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                e["acc"],
                e["type"],
                e["checked"],
                len(e["gained"]),
                len(e["lost"]),
                " | ".join(e["gained"]),
                " | ".join(e["lost"])
            ))
