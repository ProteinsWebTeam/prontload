#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import json
import logging
import math
import os
import tempfile
import time
from multiprocessing import Process, Queue

from . import oracledb


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


class MatchComparator(Process):
    def __init__(self, task_queue, result_queue, max_gap):
        super().__init__()
        self.tasks = task_queue
        self.results = result_queue
        self.max_gap = max_gap

    def run(self):
        while True:
            task = self.tasks.get()

            if task is None:
                break

            accession, matches = task
            structure = self.condense(matches, self.max_gap)
            signatures, comparisons = self.compare(matches)
            self.results.put((accession, structure, signatures, comparisons))

    @staticmethod
    def compare(matches):
        comparisons = []
        signatures = {}
        for m1 in matches:
            acc_1 = m1[0]
            len_1 = m1[2] - m1[1] + 1

            if acc_1 in signatures:
                signatures[acc_1] += 1
            else:
                signatures[acc_1] = 1

            for m2 in matches:
                acc_2 = m2[0]
                len_2 = m2[2] - m2[1] + 1

                """
                1           13      end position *is* included
                -------------       match 1 (13 - 1 + 1 = 13 aa)
                        --------    match 2 (16 - 9 + 1 = 8 aa)
                        9      16
                                    overlap = 13 - 9 + 1 = 5
                                    frac 1 = 5 / 13 = 0.38...
                                    frac 2 = 5 / 8 = 0.625
                """

                if acc_1 <= acc_2:
                    comparisons.append((acc_1, len_1, acc_2, len_2, MatchComparator.calc_overlap(m1, m2)))

        return list(signatures.items()), comparisons

    @staticmethod
    def calc_overlap(m1, m2):
        # m1 and m2 are tuples (method_ac, start, end)
        return min(m1[2], m2[2]) - max(m1[1], m2[1]) + 1

    @staticmethod
    def condense(matches, max_gap):
        # flatten matches
        locations = []

        for method_ac, pos_start, pos_end in matches:
            locations.append((pos_start, method_ac))
            locations.append((pos_end, method_ac))

        """
        Evaluate the protein's match structure, i.e. how signatures match the proteins

        -----------------------------   Protein
         <    >                         Signature 1
           <    >                       Signature 2
                      < >               Signature 3

        Flattened:
        -----------------------------   Protein
         < <  > >     < >
         1 2  1 2     3 3

        Structure, with '-' representing a "gap" (more than N bp between two positions):
        1212-33
        """

        # Sort locations by position
        locations.sort()

        """
        Do not set the offset to 0, but to the first position:
        if two proteins have the same structure,
        but the first position of one protein is > max_gap while the first position of the other protein is <= max_gap,
        a gap will be used for the first protein and not for the other, which will results in two different structures
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
    def __init__(self, task_queue, dsn, schema, **kwargs):
        super().__init__()
        self.tasks = task_queue
        self.dsn = dsn
        self.schema = schema
        self.tmpdir = kwargs.get('tmpdir', tempfile.gettempdir())
        self.chunk_size = kwargs.get('chunk_size', 100000)
        self.max_size = kwargs.get('max_size', 1000000)

        if self.tmpdir and not os.path.isdir(self.tmpdir):
            os.makedirs(self.tmpdir)

    def run(self):
        structures = {}
        signatures = {}
        comparisons = {}
        proteins = {}
        files = []

        while True:
            # Get processed proteins
            task = self.tasks.get()

            if task is None:
                break

            accession, struct, _signatures, _comparisons = task

            # struct: condensed match structure of the protein
            if struct in structures:
                code = structures[struct]
            else:
                code = structures[struct] = self.base36encode(len(structures) + 1)

            p = proteins[accession] = {
                'code': code,
                'signatures': []
            }

            # _signatures: number of times each signature matched the protein
            for acc, n_matches in _signatures:
                p['signatures'].append(acc)

                if acc not in signatures:
                    signatures[acc] = {
                        'proteins': 0,
                        'matches': 0
                    }

                signatures[acc]['proteins'] += 1
                signatures[acc]['matches'] += n_matches

            if len(proteins) == self.max_size:
                files.append(self.dump(proteins))
                proteins = {}

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

                if self.is_significant(overlap, len_1, len_2):
                    # signatures are overlapping
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
            proteins = {}

        logging.info('making predictions')

        # Get candidate signatures from DB (and non-PROSITE Pattern candidates)
        con = oracledb.connect(self.dsn)
        cur = con.cursor()
        cur.execute(
            """
            SELECT METHOD_AC, DBCODE
            FROM {}.METHOD
            WHERE CANDIDATE = 'Y'
            """.format(self.schema)
        )

        candidates = set()
        non_prosite_candidates = set()
        for accession, dbcode in cur:
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

                    adj_1 = adj_relations.get(acc_2, {}).get(acc_1,
                                                             0)  # inverting acc2 and acc1 to have the same predictions than HH, todo: check if it is really ok
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

        # Inserting predictions and signature stats
        logging.info('creating METHOD_MATCH table')
        oracledb.drop_table(cur, self.schema, 'METHOD_MATCH')
        cur.execute(
            """
            CREATE TABLE {}.METHOD_MATCH
            (
                METHOD_AC VARCHAR2(25) NOT NULL,
                N_MATCHES NUMBER NOT NULL,
                N_PROT NUMBER NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        # Populating METHOD_MATCH
        signatures = [(acc, s['matches'], s['proteins']) for acc, s in signatures.items()]
        for i in range(0, len(signatures), self.chunk_size):
            cur.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.METHOD_MATCH (METHOD_AC, N_MATCHES, N_PROT)
                VALUES (:1, :2, :3)
                """.format(self.schema),
                signatures[i:i + self.chunk_size]
            )
            con.commit()

        # Primary key for METHOD_MATCH
        cur.execute(
            """
            ALTER TABLE {}.METHOD_MATCH
            ADD CONSTRAINT PK_METHOD_MATCH PRIMARY KEY (METHOD_AC)
            """.format(self.schema)
        )

        # Stats for METHOD_MATCH
        oracledb.gather_stats(cur, self.schema, 'METHOD_MATCH', cascade=True)

        # Privileges
        oracledb.grant(cur, 'SELECT', self.schema, 'METHOD_MATCH', 'INTERPRO_SELECT')

        logging.info('creating METHOD_OVERLAP table')
        oracledb.drop_table(cur, self.schema, 'METHOD_OVERLAP')
        cur.execute(
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
            cur.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.METHOD_OVERLAP (
                    METHOD_AC1, METHOD_AC2, N_PROT, N_OVER, N_PROT_OVER, AVG_OVER, AVG_FRAC1, AVG_FRAC2
                )
                VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
                """.format(self.schema),
                overlaps[i:i + self.chunk_size]
            )
            con.commit()

        # Primary key for METHOD_OVERLAP
        cur.execute(
            """
            ALTER TABLE {}.METHOD_OVERLAP
            ADD CONSTRAINT PK_METHOD_OVERLAP PRIMARY KEY (METHOD_AC1, METHOD_AC2)
            """.format(self.schema)
        )

        # Stats for METHOD_OVERLAP
        oracledb.gather_stats(cur, self.schema, 'METHOD_OVERLAP', cascade=True)

        # Privileges
        oracledb.grant(cur, 'SELECT', self.schema, 'METHOD_OVERLAP', 'INTERPRO_SELECT')

        logging.info('creating METHOD_PREDICTION table')
        oracledb.drop_table(cur, self.schema, 'METHOD_PREDICTION')
        cur.execute(
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
            cur.executemany(
                """
                INSERT /*+APPEND*/ INTO {}.METHOD_PREDICTION (METHOD_AC1, METHOD_AC2, RELATION)
                VALUES (:1, :2, :3)
                """.format(self.schema),
                predictions[i:i + self.chunk_size]
            )
            con.commit()

        # Primary key for METHOD_PREDICTION
        cur.execute(
            """
            ALTER TABLE {}.METHOD_PREDICTION
            ADD CONSTRAINT PK_METHOD_PREDICTION PRIMARY KEY (METHOD_AC1, METHOD_AC2)
            """.format(self.schema)
        )

        # Stats for METHOD_PREDICTION
        oracledb.gather_stats(cur, self.schema, 'METHOD_PREDICTION', cascade=True)

        # Privileges
        oracledb.grant(cur, 'SELECT', self.schema, 'METHOD_PREDICTION', 'INTERPRO_SELECT')

        # Staging table (store the condensation code for proteins, stored in files at the moment)
        logging.info('creating METHOD2PROTEIN_STG table')
        oracledb.drop_table(cur, self.schema, 'METHOD2PROTEIN_STG')
        cur.execute(
            """
            CREATE TABLE {}.METHOD2PROTEIN_STG
            (
                METHOD_AC VARCHAR2(25) NOT NULL,
                PROTEIN_AC VARCHAR2(15) NOT NULL,
                CONDENSE VARCHAR(100) NOT NULL
            ) NOLOGGING
            """.format(self.schema)
        )

        # Populating METHOD2PROTEIN_STG by loading files
        for filepath in files:
            with gzip.open(filepath, 'rt') as fh:
                proteins = json.load(fh)

            os.unlink(filepath)

            data = [(m_ac, p_ac, p['code']) for p_ac, p in proteins.items() for m_ac in p['signatures']]
            for i in range(0, len(data), self.chunk_size):
                cur.executemany(
                    """
                    INSERT /*+APPEND*/ INTO {}.METHOD2PROTEIN_STG (METHOD_AC, PROTEIN_AC, CONDENSE)
                    VALUES (:1, :2, :3)
                    """.format(self.schema),
                    data[i:i + self.chunk_size]
                )
                con.commit()

            # Can I haz some memory back?
            proteins = None
            data = None

        # Indexing METHOD2PROTEIN_STG (no primary key)
        cur.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN_STG$PROTEIN
            ON {}.METHOD2PROTEIN_STG (PROTEIN_AC)
            NOLOGGING
            """.format(self.schema)
        )

        # Stats for METHOD2PROTEIN_STG
        oracledb.gather_stats(cur, self.schema, 'METHOD2PROTEIN_STG', cascade=True)

        logging.info('creating METHOD2PROTEIN table')
        oracledb.drop_table(cur, self.schema, 'METHOD2PROTEIN')
        cur.execute(
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

        # Populating METHOD2PROTEIN (merging METHOD2PROTEIN_STG, PROTEIN, AND ETAXI)
        cur.execute(
            """
            INSERT /*+APPEND*/ INTO {0}.METHOD2PROTEIN (METHOD_AC, PROTEIN_AC, DBCODE, CONDENSE, LEN, LEFT_NUMBER, DESC_ID)
            SELECT M2P.METHOD_AC, M2P.PROTEIN_AC, P.DBCODE, M2P.CONDENSE, P.LEN, NVL(E.LEFT_NUMBER, 0), PD.DESC_ID
            FROM {0}.METHOD2PROTEIN_STG M2P
            INNER JOIN {0}.PROTEIN P ON M2P.PROTEIN_AC = P.PROTEIN_AC
            INNER JOIN {0}.ETAXI E ON P.TAX_ID = E.TAX_ID
            INNER JOIN {0}.PROTEIN_DESC PD ON M2P.PROTEIN_AC = PD.PROTEIN_AC
            """.format(self.schema)
        )
        con.commit()

        # Dropping staging table
        oracledb.drop_table(cur, self.schema, 'METHOD2PROTEIN_STG')

        # Indexing METHOD2PROTEIN
        cur.execute(
            """
            ALTER TABLE {}.METHOD2PROTEIN
            ADD CONSTRAINT PK_METHOD2PROTEIN PRIMARY KEY (METHOD_AC, PROTEIN_AC)
            """.format(self.schema)
        )
        cur.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$M
            ON {}.METHOD2PROTEIN (METHOD_AC)
            NOLOGGING
            """.format(self.schema)
        )
        cur.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$P
            ON {}.METHOD2PROTEIN (PROTEIN_AC)
            NOLOGGING
            """.format(self.schema)
        )
        cur.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$LN
            ON {}.METHOD2PROTEIN (LEFT_NUMBER)
            NOLOGGING
            """.format(self.schema)
        )
        cur.execute(
            """
            CREATE INDEX I_METHOD2PROTEIN$M$DB
            ON {}.METHOD2PROTEIN (METHOD_AC, DBCODE)
            NOLOGGING
            """.format(self.schema)
        )

        # Stats for METHOD2PROTEIN
        oracledb.gather_stats(cur, self.schema, 'METHOD2PROTEIN', cascade=True)

        # Privileges
        oracledb.grant(cur, 'SELECT', self.schema, 'METHOD2PROTEIN', 'INTERPRO_SELECT')

        cur.close()
        con.close()

    def dump(self, proteins):
        fd, filepath = tempfile.mkstemp(suffix='.json.gz', dir=self.tmpdir)
        os.close(fd)

        with gzip.open(filepath, 'wt') as fh:
            json.dump(proteins, fh)

        return filepath

    @staticmethod
    def base36encode(integer):
        chars = '0123456789abcdefghijklmnopqrstuvwxyz'
        encoded = ''

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            encoded = chars[remainder] + encoded

        return encoded

    @staticmethod
    def is_significant(o, l1, l2):
        # Consider that matches significantly overlap if the overlap is longer than the half of the shortest match
        return o > (min(l1, l2) / 2)


def load_matches_and_predictions(dsn, schema, **kwargs):
    max_gap = kwargs.get('max_gap', 20)
    processes = kwargs.get('processes', 3)
    chunk_size = kwargs.get('chunk_size', 100000)
    max_size = kwargs.get('max_size', 10000000)

    proteins = Queue(maxsize=100000)
    comparisons = Queue(maxsize=100000)

    comparators = [
        MatchComparator(proteins, comparisons, max_gap)
        for _ in range(max(1, processes-2))  # minus two processes: parent process + MatchAggregator
    ]

    aggregator = MatchAggregator(comparisons, dsn, schema, **kwargs)

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
    logging.info('dropping / creating MATCH table')
    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'MATCH')
    cur.execute(
        """
        CREATE TABLE {}.MATCH
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            MODEL_AC VARCHAR2(25) NOT NULL,
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            STATUS CHAR(1) NOT NULL,
            DBCODE CHAR(1) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    logging.info('processing InterPro matches')
    cur.execute(
        """
        SELECT
          MA.PROTEIN_AC,
          MA.METHOD_AC,
          MA.MODEL_AC,
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
        model_ac = row[1] if row[2] is None else row[2]
        pos_start = math.floor(row[3])
        pos_end = math.floor(row[4])
        status = row[5]
        dbcode = row[6]
        method_type = row[7]
        protein_length = math.floor(row[8])
        is_fragment = row[9] == 'Y'

        matches_all.append((protein, method_ac, model_ac, pos_start, pos_end, status, dbcode))

        if len(matches_all) == max_size:
            cnt_matches += max_size
            for i in range(0, len(matches_all), chunk_size):
                cur2.executemany(
                    """
                    INSERT /*+APPEND*/ INTO {}.MATCH (PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, STATUS, DBCODE)
                    VALUES (:1, :2, :3, :4, :5, :6, :7)
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
            INSERT /*+APPEND*/ INTO {}.MATCH (PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, STATUS, DBCODE)
            VALUES (:1, :2, :3, :4, :5, :6, :7)
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

    # Get descriptions for before update (METHOD2SWISS_DE was populated during protein update)
    con = oracledb.connect(dsn)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, M.DESCRIPTION
        FROM {0}.METHOD2SWISS_DE M 
        INNER JOIN {0}.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        """.format(schema)
    )

    old_entries = {}
    for acc, descr in cur:
        if acc in old_entries:
            old_entries[acc].append(descr)
        else:
            old_entries[acc] = [descr]

    # Get descriptions for after update
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
                'gained': [descr for descr in new_entry_descrs if descr not in old_entries[acc]]
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
