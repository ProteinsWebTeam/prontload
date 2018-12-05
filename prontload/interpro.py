#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bisect
import hashlib
import heapq
import logging
import os
import time
from multiprocessing import Process, Queue
from typing import Tuple

from . import io
from .oracledb import Connection


BULK_INSERT_SIZE = 100000


def create_synonyms(dsn, src, dst):
    tables = (
        "ENTRY",
        "ENTRY2METHOD",
        "ENTRY2ENTRY",
        "ENTRY2COMP",
        "METHOD2SWISS_DE",
        "METHOD_COMMENT",
        "USER_PRONTO"
    )

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
        INSERT /*+APPEND*/ INTO {}.CV_DATABASE (
            DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE
        )
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        """.format(schema)
    )
    con.commit()

    con.optimise_table(schema, "CV_DATABASE", cascade=True)
    con.grant("SELECT", schema, "CV_DATABASE", "INTERPRO_SELECT")


def load_matches(dsn, schema):
    con = Connection(dsn)

    # Dropping (if exists) and recreating table
    con.drop_table(schema, "MATCH")
    con.execute(
        """
        CREATE TABLE {}.MATCH
        (
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            METHOD_AC VARCHAR2(25) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            MODEL_AC VARCHAR2(25),
            POS_FROM NUMBER(5) NOT NULL,
            POS_TO NUMBER(5) NOT NULL,
            FRAGMENTS VARCHAR2(200) DEFAULT NULL
        ) NOLOGGING
        """.format(schema)
    )

    """
    Insert directly signature matches and MobiDB-lite matches
    For some signature matches, the HMM model accession
        is the signature accession itself: in such cases, use NULL
    """
    con.execute(
        """
        INSERT /*+APPEND*/ INTO {}.MATCH
        SELECT
          PROTEIN_AC,
          METHOD_AC,
          DBCODE,
          CASE
            WHEN MODEL_AC IS NOT NULL AND MODEL_AC != METHOD_AC
            THEN MODEL_AC
            ELSE NULL
          END AS MODEL_AC,
          POS_FROM,
          POS_TO,
          FRAGMENTS
        FROM INTERPRO.MATCH
        UNION ALL
        SELECT
          PROTEIN_AC,
          METHOD_AC,
          DBCODE,
          NULL,
          POS_FROM,
          POS_TO,
          NULL
        FROM INTERPRO.FEATURE_MATCH
        WHERE DBCODE = 'g'
        """.format(schema)
    )
    con.commit()

    # Finalizing table
    con.execute(
        """
        CREATE INDEX I_MATCH$PROTEIN
        ON {}.MATCH (PROTEIN_AC) NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_MATCH$METHOD
        ON {}.MATCH (METHOD_AC) NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_MATCH$DBCODE
        ON {}.MATCH (DBCODE) NOLOGGING
        """.format(schema)
    )
    con.optimise_table(schema, "MATCH", cascade=True)
    con.grant("SELECT", schema, "MATCH", "INTERPRO_SELECT")

    # TODO: ensure that `METHOD` is ready
    con.execute(
        """
        UPDATE {0}.METHOD ME
        SET PROTEIN_COUNT = (
            SELECT COUNT(DISTINCT PROTEIN_AC)
            FROM {0}.MATCH MA
            WHERE ME.METHOD_AC = MA.METHOD_AC
        )
        """.format(schema)
    )
    con.commit()
    cur.close()


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

    con.optimise_table(schema, "PROTEIN", cascade=True)
    con.grant("SELECT", schema, "PROTEIN", "INTERPRO_SELECT")


def load_signatures(dsn, schema):
    con = Connection(dsn)
    con.drop_table(schema, "METHOD")
    con.execute(
        """
        CREATE TABLE {}.METHOD
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            CANDIDATE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(220),
            SIG_TYPE CHAR(1),
            PROTEIN_COUNT NUMBER(8) DEFAULT 0 NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    con.execute(
        """
        INSERT /*+APPEND*/ INTO {}.METHOD (
            METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE
        )
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

    con.optimise_table(schema, "METHOD", cascade=True)
    con.grant("SELECT", schema, "METHOD", "INTERPRO_SELECT")


def load_taxa(dsn, schema):
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

    con.optimise_table(schema, "ETAXI", cascade=True)
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

    for i in range(0, len(lineage), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.LINEAGE (LEFT_NUMBER, TAX_ID, RANK)
            VALUES (:1, :2, :3)
            """.format(schema),
            lineage[i:i+BULK_INSERT_SIZE]
        )
        con.commit()

    con.execute(
        """
        CREATE INDEX I_LINEAGE$L$R
        ON {}.LINEAGE (LEFT_NUMBER, RANK)
        NOLOGGING
        """.format(schema)
    )

    con.optimise_table(schema, "LINEAGE", cascade=True)
    con.grant("SELECT", schema, "LINEAGE", "INTERPRO_SELECT")


def chunk_signatures(con, schema, size=1000):
    accessions = []
    cnt = size
    for row in con.get(
            """
                SELECT METHOD_AC
                FROM {}.METHOD
                ORDER BY METHOD_AC
            """.format(schema)
    ):
        if cnt == size:
            accessions.append(row[0])
            cnt = 1
        else:
            cnt += 1

    return accessions


class ProteinConsumer(Process):
    def __init__(self, dsn, schema, max_gap, queue_in, queue_out, dir=None):
        super().__init__()
        self.dsn = dsn
        self.schema = schema
        self.max_gap = max_gap
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.dir = dir

    def run(self):
        con = Connection(self.dsn)

        # Get signatures
        accessions = chunk_signatures(con, self.schema)
        organiser_names = io.Organiser(accessions, dir=self.dir)
        organiser_taxa = io.Organiser(accessions, dir=self.dir)

        # Get lineages for the METHOD_TAXA table
        ranks = {
            "superkingdom", "kingdom", "phylum", "class", "order",
            "family", "genus", "species"
        }
        left_numbers = {}
        for rank, left_num, tax_id in con.get(
                """
                SELECT RANK, LEFT_NUMBER, TAX_ID
                FROM {}.LINEAGE
                """.format(self.schema)
        ):
            if rank not in ranks:
                continue
            elif left_num not in left_numbers:
                left_numbers[left_num] = {}
            left_numbers[left_num][rank] = tax_id

        signatures = {}
        comparisons = {}
        residue_coverages = {}
        residue_overlaps = {}
        while True:
            chunk = self.queue_in.get()
            if chunk is None:
                break

            data = []
            for args in chunk:
                protein_acc, dbcode, length, desc_id, left_num, matches = args

                md5 = self.hash(matches, self.max_gap)
                _signatures, _comparisons = self.compare(matches)
                ranks = left_numbers.get(left_num, {"no rank": -1})

                for method_acc in _signatures:
                    data.append((
                        method_acc, protein_acc, dbcode,
                        md5, length, left_num, desc_id
                    ))

                    if method_acc not in signatures:
                        signatures[method_acc] = {
                            "proteins": 0,
                            "matches": 0
                        }

                    signatures[method_acc]["proteins"] += 1
                    n_matches = _signatures[method_acc]
                    signatures[method_acc]["matches"] += n_matches

                    # Taxonomic origins
                    for rank, tax_id in ranks.items():
                        organiser_taxa.add(method_acc, (rank, tax_id))

                    # UniProt descriptions
                    organiser_names.add(method_acc, (desc_id, dbcode))

                # _comparisons: match overlaps between signatures
                for acc_1 in _comparisons:
                    if acc_1 in comparisons:
                        d = comparisons[acc_1]
                    else:
                        d = comparisons[acc_1] = {}

                    for acc_2 in _comparisons[acc_1]:
                        if acc_2 in d:
                            comp = d[acc_2]
                        else:
                            comp = d[acc_2] = {
                                # num of proteins in which both sign. occur
                                # (not necessarily overlapping)
                                "prot": 0,
                                # num of proteins in which signatures overlap
                                "prot_over": 0,
                                # num of times signatures overlap
                                # (>= prot_over)
                                "over": 0,
                                # sum of overlap lengths
                                # (to compute the average length later)
                                "length": 0,
                                # sum of fractions of matches overlapping
                                # (overlap length / match length)
                                "frac_1": 0,
                                "frac_2": 0
                            }

                        comp["prot"] += 1
                        prot_over = False
                        for l1, l2, o in _comparisons[acc_1][acc_2]:
                            if o > (min(l1, l2) / 2):
                                """
                                Consider that matches significantly overlap
                                if the overlap is longer
                                than the half of the shortest match
                                """
                                comp["frac_1"] += o / l1
                                comp["frac_2"] += o / l2
                                comp["over"] += 1
                                comp["length"] += o
                                prot_over = True

                        if prot_over:
                            comp["prot_over"] += 1

                self.compare_fragments(matches, residue_coverages,
                                       residue_overlaps)

            for i in range(0, len(data), BULK_INSERT_SIZE):
                con.executemany(
                    """
                    INSERT /*+APPEND*/ INTO {}.METHOD2PROTEIN (
                        METHOD_AC, PROTEIN_AC, DBCODE,
                        MD5, LEN, LEFT_NUMBER, DESC_ID
                    )
                    VALUES (:1, :2, :3, :4, :5, :6, :7)
                    """.format(self.schema),
                    data[i:i+BULK_INSERT_SIZE]
                )
                con.commit()

            organiser_names.dump()
            organiser_taxa.dump()

        size = organiser_names.merge()
        logging.info("names: {} bytes".format(size))

        size = organiser_taxa.merge()
        logging.info("taxa: {} bytes".format(size))

        self.queue_out.put((
            signatures,
            comparisons,
            residue_coverages,
            residue_overlaps,
            organiser_names,
            organiser_taxa,
        ))

    @staticmethod
    def compare(matches):
        comparisons = {}
        signatures = {}

        # Fourth item is fragments_str, which is not used here
        for m1_acc, m1_start, m1_end, _ in matches:
            m1_len = m1_end - m1_start + 1

            if m1_acc in signatures:
                signatures[m1_acc] += 1
                m1 = comparisons[m1_acc]
            else:
                signatures[m1_acc] = 1
                m1 = comparisons[m1_acc] = {}

            for m2_acc, m2_start, m2_end, _ in matches:
                if m1_acc > m2_acc:
                    continue
                elif m2_acc in m1:
                    m2 = m1[m2_acc]
                else:
                    m2 = m1[m2_acc] = []

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

                # overlap
                m2.append((
                    m1_len,
                    m2_len,
                    min(m1_end, m2_end) - max(m1_start, m2_start) + 1
                ))

        return signatures, comparisons

    @staticmethod
    def hash(matches, max_gap):
        # flatten matches
        locations = []

        # Fourth item is fragments_str, which is not used here
        for method_ac, pos_start, pos_end, _ in matches:
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

        # overall match structure
        structure = []
        # close signatures (less than max_gap between two positions)
        methods = []

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

        return hashlib.md5(
            '/'.join(structure).encode("utf-8")
        ).hexdigest()

    @staticmethod
    def compare_fragments(matches, residues, overlaps):
        # Parse the fragments string
        fragments = []
        for method_acc, start, end, fragments_str in matches:
            if fragments_str is None:
                fragments.append((method_acc, start, end))
            else:
                # Discontinuous domains
                fallback = True
                for frag in fragments_str.split(','):
                    """
                    Format: START-END-TYPE
                    Types:
                        * S: Continuous single chain domain
                        * N: N-terminal discontinuous
                        * C: C-terminal discontinuous
                        * NC: N and C -terminal discontinuous
                    """
                    pos_start, pos_end, t = frag.split('-')
                    pos_start = int(pos_start)
                    pos_end = int(pos_end)
                    if pos_start < pos_end:
                        # At least one well formated discontinuous domain
                        fallback = False
                        fragments.append((method_acc, pos_start, pos_end))

                if fallback:
                    # Fallback to match start/end positions
                    fragments.append((method_acc, start, end))

        # Sort by positions
        fragments.sort(key=lambda x: (x[1], x[2]))

        # Group sorted fragments by signature
        signatures = {}
        for method_acc, start, end in fragments:
            if method_acc in signatures:
                signatures[method_acc].append((start, end))
            else:
                signatures[method_acc] = [(start, end)]

        for acc_1 in signatures:
            # number of residues covered by the signature matches
            _residues = sum([e - s + 1 for s, e in signatures[acc_1]])

            if acc_1 in residues:
                residues[acc_1] += _residues
            else:
                residues[acc_1] = _residues
                overlaps[acc_1] = {}

            for acc_2 in signatures:
                if acc_1 < acc_2:
                    if acc_2 not in overlaps[acc_1]:
                        overlaps[acc_1][acc_2] = 0

                    i = 0
                    start_2, end_2 = signatures[acc_2][i]
                    for start_1, end_1 in signatures[acc_1]:
                        while end_2 < start_1:
                            i += 1
                            try:
                                start_2, end_2 = signatures[acc_2][i]
                            except IndexError:
                                break

                        o = min(end_1, end_2) - max(start_1, start_2) + 1
                        if o > 0:
                            overlaps[acc_1][acc_2] += o


def dump_proteins(dsn, schema, dst):
    with io.Store(dst) as store:
        con = Connection(dsn)
        cnt = 0
        for row in con.get(
                """
                SELECT
                  P.PROTEIN_AC, P.LEN, P.DBCODE, D.DESC_ID,
                  NVL(E.LEFT_NUMBER, 0)
                FROM {0}.PROTEIN P
                INNER JOIN {0}.ETAXI E
                  ON P.TAX_ID = E.TAX_ID
                INNER JOIN {0}.PROTEIN_DESC D
                  ON P.PROTEIN_AC = D.PROTEIN_AC
                WHERE P.FRAGMENT = 'N'
                ORDER BY P.PROTEIN_AC
                """.format(schema)
        ):
            store.add(row)
            cnt += 1
            if not cnt % 10000000:
                logging.info("proteins: {:>12}".format(cnt))

        logging.info("proteins: {:>12}".format(cnt))


def dump_matches(dsn, schema, processes, dst, dir=None, bucket_size=1000000,
                 buffer_size=1000000):
    con = Connection(dsn)
    i = bucket_size
    keys = []
    for accession, in con.get(
            """
            SELECT PROTEIN_AC
            FROM {}.PROTEIN
            WHERE FRAGMENT = 'N'
            ORDER BY PROTEIN_AC
            """.format(schema)
    ):
        if i == bucket_size:
            keys.append(accession)
            i = 1
        else:
            i += 1

    out_queue = Queue()
    processes = max(1, processes - 1)
    organisers = []
    queues = []
    workers = []
    chunks = []
    _keys = []
    step = int(len(keys) / processes + 1)
    for i in range(0, len(keys), step):
        o = io.Organiser(keys[i:i+step], dir=dir)
        q = Queue(maxsize=1)
        w = Process(target=organise_matches, args=(o, q, out_queue))
        w.start()
        organisers.append(o)
        queues.append(q)
        workers.append(w)
        chunks.append([])
        _keys.append(keys[i])

    cnt = 0
    keys = _keys
    buffer_size //= processes
    for row in con.get(
            """
            SELECT
              MA.PROTEIN_AC,
              MA.METHOD_AC,
              MA.DBCODE,
              ME.SIG_TYPE,
              MA.POS_FROM,
              MA.POS_TO,
              MA.FRAGMENTS
            FROM INTERPRO.MATCH MA
              INNER JOIN {0}.METHOD ME
                ON MA.METHOD_AC = ME.METHOD_AC
            WHERE MA.PROTEIN_AC IN (
              SELECT PROTEIN_AC
              FROM {0}.PROTEIN
              WHERE FRAGMENT = 'N'
            )
            """.format(schema)
    ):
        protein_acc = row[0]
        method_acc = row[1]
        method_dbcode = row[2]
        method_type = row[3]
        pos_start = int(row[4])
        pos_end = int(row[5])
        fragments_str = row[6]

        # Find worker feeding the organiser of this protein
        i = bisect.bisect_right(keys, protein_acc) - 1
        chunks[i].append((protein_acc, method_acc, method_dbcode, method_type,
                          pos_start, pos_end, fragments_str))

        if len(chunks[i]) == buffer_size:
            queues[i].put(chunks[i])
            chunks[i] = []

        cnt += 1
        if not cnt % 100000000:
            logging.info("matches: {:>12}".format(cnt))

    # Trigger organiser merging
    for c, q in zip(chunks, queues):
        q.put(c)
        q.put(None)

    logging.info("matches: {:>12}".format(cnt))

    # Get temporary space used by each organiser
    size = sum([out_queue.get() for _ in workers])

    for w in workers:
        w.join()

    logging.info("organisers disk space: {} bytes".format(size))

    with io.Store(dst) as store:
        for o in organisers:
            for protein_acc, matches in o:
                store.add((protein_acc, matches))

            o.remove()

        logging.info("store disk space: {} bytes".format(store.size))


def organise_matches(organiser: io.Organiser, in_queue: Queue,
                     out_queue: Queue):
    while True:
        chunk = in_queue.get()
        if chunk is None:
            break

        for accession, *match in chunk:
            organiser.add(accession, match)

        organiser.dump()

    size = organiser.merge()
    out_queue.put(size)


def _process_proteins(con, dsn, schema, organisers, store, max_gap,
                     chunk_size=10000, tmpdir=None):
    queue_in = Queue(maxsize=len(organisers))
    queue_out = Queue()

    workers = []
    for _ in organisers:
        w = ProteinConsumer(dsn, schema, max_gap, queue_in, queue_out,
                            dir=tmpdir)
        w.start()
        workers.append(w)

    cnt = 0
    chunk = []
    ts = time.time()

    # Open store for reading
    iter(store)

    # Iterate by organiser
    for organiser in organisers:
        # Iterate by the organiser's items
        for protein_acc, matches in organiser:
            # Read store until we find the current protein
            acc, length, prot_dbcode, desc_id, left_num = next(store)
            while acc < protein_acc:
                acc, length, prot_dbcode, desc_id, left_num = next(store)

            if acc != protein_acc:
                # Protein missing from store
                raise KeyError("{}: missing protein".format(protein_acc))

            _matches = []
            methods = {}
            methods_fragments = {}
            for m in matches:
                (method_acc, method_dbcode, method_type,
                 pos_start, pos_end, fragments_str) = m

                if method_dbcode in ('F', 'V'):
                    """
                    PANTHER & PRINTS:
                        Merge protein matches.
                        If the signature is a family*,
                            use the entire protein.

                        * all PANTHER signatures,
                            almost all PRINTS signatures
                    """
                    if method_acc not in methods:
                        methods_fragments[method_acc] = fragments_str

                        if method_type == 'F':
                            # Families: use the entire protein sequence
                            methods[method_acc] = [(1, length)]
                        else:
                            methods[method_acc] = [(pos_start, pos_end)]
                    elif method_type != 'F':
                        # Only for non-families
                        # (because families use the entire protein)
                        methods[method_acc].append((pos_start, pos_end))
                else:
                    _matches.append((method_acc, pos_start, pos_end,
                                     fragments_str))

            # Merge matches
            for method_acc in methods:
                min_pos = None
                max_pos = None

                for pos_start, pos_end in methods[method_acc]:
                    if min_pos is None or pos_start < min_pos:
                        min_pos = pos_start
                    if max_pos is None or pos_end > max_pos:
                        max_pos = pos_end

                _matches.append((method_acc, min_pos, max_pos,
                                 methods_fragments[method_acc]))

            chunk.append((protein_acc, prot_dbcode, length, desc_id,
                          left_num, _matches))

            if len(chunk) == chunk_size:
                # Enqueue chunk of proteins for ProteinConsumers
                queue_in.put(chunk)
                chunk = []

            cnt += 1
            if not cnt % 1000000:
                logging.info("{:>12} ({:.0f} proteins/sec)".format(
                    cnt,
                    cnt / (time.time() - ts)
                ))

    if chunk:
        queue_in.put(chunk)
        chunk = []

    for _ in workers:
        queue_in.put(None)

    logging.info("{:>12} ({:.0f} proteins/sec)".format(
        cnt,
        cnt / (time.time() - ts)
    ))

    # For predictions
    signatures = {}
    comparisons = {}
    # For residue overlap comparison
    residue_coverages = {}
    residue_overlaps = {}
    # For counts
    name_organisers = []
    taxon_organisers = []
    for _ in workers:
        s, c, rc, ro, no, to = queue_out.get()

        # Merge dictionaries
        for acc in s:
            if acc in signatures:
                for k, v in s[acc].items():
                    signatures[acc][k] += v
            else:
                signatures[acc] = s[acc]

        for acc_1 in c:
            if acc_1 in comparisons:
                d = comparisons[acc_1]
            else:
                comparisons[acc_1] = c[acc_1]
                continue

            for acc_2 in c[acc_1]:
                if acc_2 in d:
                    for k, v in c[acc_1][acc_2].items():
                        d[acc_2][k] += v
                else:
                    d[acc_2] = c[acc_1][acc_2]

        for acc in rc:
            if acc in residue_coverages:
                residue_coverages[acc] += rc[acc]
            else:
                residue_coverages[acc] = rc[acc]

        for acc_1 in ro:
            if acc_1 in residue_overlaps:
                for acc_2 in ro[acc_1]:
                    if acc_2 in residue_overlaps[acc_1]:
                        residue_overlaps[acc_1][acc_2] += ro[acc_1][acc_2]
                    else:
                        residue_overlaps[acc_1][acc_2] = ro[acc_1][acc_2]
            else:
                residue_overlaps[acc_1] = ro[acc_1]

        name_organisers.append(no)
        taxon_organisers.append(to)

    for w in workers:
        w.join()

    con.drop_table(schema, "METHOD_SIMILARITY")
    con.execute(
        """
        CREATE TABLE {}.METHOD_SIMILARITY
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            SIMILARITY BINARY_DOUBLE NOT NULL,
            CONTAINMENT1 BINARY_DOUBLE NOT NULL,
            CONTAINMENT2 BINARY_DOUBLE NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    data = []
    for acc_1 in residue_overlaps:
        for acc_2 in residue_overlaps[acc_1]:
            i = residue_overlaps[acc_1][acc_2]
            u = residue_coverages[acc_1] + residue_coverages[acc_2] - i
            index = i / u
            cont_1 = i / residue_coverages[acc_1]
            cont_2 = i / residue_coverages[acc_2]
            data.append((acc_1, acc_2, index, cont_1, cont_2))

    for i in range(0, len(data), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT /*+APPEND*/
            INTO {}.METHOD_SIMILARITY
            VALUES (:1, :2, :3, :4, :5)
            """.format(schema),
            data[i:i+BULK_INSERT_SIZE]
        )
        con.commit()

    con.execute(
        """
        ALTER TABLE {}.METHOD_SIMILARITY
        ADD CONSTRAINT PK_METHOD_SIMILARITY
        PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        """.format(schema)
    )

    return signatures, comparisons, name_organisers, taxon_organisers


def make_predictions(dsn, schema, signatures, comparisons):
    candidates = set()
    non_prosite_candidates = set()
    con = Connection(dsn)
    for method_acc, dbcode in con.get(
            """
            SELECT METHOD_AC, DBCODE
            FROM {}.METHOD
            WHERE CANDIDATE = 'Y'
            """.format(schema)
    ):
        candidates.add(method_acc)
        if dbcode != 'P':
            non_prosite_candidates.add(method_acc)

    # Determine relations/adjacent
    relations = {}
    adjacents = {}
    for acc_1 in comparisons:
        s_1 = signatures[acc_1]

        for acc_2 in comparisons[acc_1]:
            s_2 = signatures[acc_2]
            c = comparisons[acc_1][acc_2]

            if (c["prot_over"] >= s_1["proteins"] * 0.4
                    or c["prot_over"] >= s_2["proteins"] * 0.4):
                # is a relation
                if acc_1 in relations:
                    relations[acc_1].append(acc_2)
                else:
                    relations[acc_1] = [acc_2]
            elif c["prot"] >= s_1["proteins"] * 0.1:
                # is adjacent
                if acc_1 in adjacents:
                    adjacents[acc_1].append(acc_2)
                else:
                    adjacents[acc_1] = [acc_2]

            if acc_1 != acc_2:
                if (c["prot_over"] >= s_1["proteins"] * 0.4
                        or c["prot_over"] >= s_2["proteins"] * 0.4):
                    # is a relation
                    if acc_2 in relations:
                        relations[acc_2].append(acc_1)
                    else:
                        relations[acc_2] = [acc_1]
                elif c["prot"] >= s_2["proteins"] * 0.1:
                    # is adjacent
                    if acc_2 in adjacents:
                        adjacents[acc_2].append(acc_1)
                    else:
                        adjacents[acc_2] = [acc_1]

    """
    Determine extra/adjacent relation:

    Let A, and B be two signatures in a relationship.
    extra relation:
        num of signatures in a relationship with A but not with B
    adjacent relation:
        num of signatures in a relationship with A and adjacent to B
    """
    extra_relations = {}
    adj_relations = {}
    for acc_1 in relations:
        # signatures in a relationship with acc_1
        # (that are candidate and not from PROSITE pattern)
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

            if (c["prot_over"] >= s_1["proteins"] * 0.4
                    or c["prot_over"] >= s_2["proteins"] * 0.4):
                extra_1 = extra_relations.get(acc_1, {}).get(acc_2, 0)
                extra_2 = extra_relations.get(acc_2, {}).get(acc_1, 0)

                """
                Inverting acc_2 and acc_1
                to have the same predictions than HH
                """
                # TODO: is this really OK?
                adj_1 = adj_relations.get(acc_2, {}).get(acc_1, 0)
                adj_2 = adj_relations.get(acc_1, {}).get(acc_2, 0)

                if c["over"]:
                    """
                    frac_1 and frac_2 are the sum of ratios
                    (overlap / match length).
                    if close to 1:
                        the match was mostly contained
                        by the overlap in most cases

                            A   -------------
                            B       ---

                        frac(B) = overlap(A,B) / length(B)
                                = 1
                                indeed B is 100% within the overlap

                    len_1 and len_2 are the average of sums.
                    If len(B) is close to 1:
                        B was mostly within the overlap in most cases
                        so B < A (because B ~ overlap and A >= overlap)
                        so B CONTAINED_BY A
                    """
                    len_1 = c["frac_1"] / c["over"]
                    len_2 = c["frac_2"] / c["over"]
                else:
                    len_1 = len_2 = 0

                """
                Parent/Child relationships:
                The protein/matches made by the child entry
                    must be a complete (>75%) subset of the parent entry

                if over(A) > 0.75, it means that A overlaps with B
                    in at least 75% of its proteins or matches:
                        A is a CHILD_OF B
                """
                over_1 = min(
                    c["over"] / s_1["matches"],
                    c["prot_over"] / s_1["proteins"]
                )
                over_2 = min(
                    c["over"] / s_2["matches"],
                    c["prot_over"] / s_2["proteins"]
                )
                if len_1 >= 0.5 and len_2 >= 0.5:
                    if (over_1 > 0.75 and over_2 >= 0.75 and
                            not extra_1 and not extra_2):
                        prediction = "ADD_TO"
                    elif over_1 > 0.75 and not extra_1 and not adj_1:
                        prediction = "CHILD_OF"
                    elif over_2 > 0.75 and not extra_2 and not adj_2:
                        prediction = "PARENT_OF"  # acc_2 child of acc_1
                    elif len_1 >= 0.9:
                        if len_2 >= 0.9:
                            prediction = "C/C"
                        else:
                            prediction = "CONTAINED_BY"
                    elif len_2 >= 0.9:
                        prediction = "CONTAINER_OF"
                    else:
                        prediction = "OVERLAPS"
                elif len_1 >= 0.9:
                    prediction = "CONTAINED_BY"
                elif len_2 >= 0.9:
                    prediction = "CONTAINER_OF"
                else:
                    prediction = "OVERLAPS"

                predictions.append((acc_1, acc_2, prediction))

                # switch (acc_1, acc_2) -> (acc_2, acc_1)
                if prediction == "CHILD_OF":
                    prediction = "PARENT_OF"
                elif prediction == "PARENT_OF":
                    prediction = "CHILD_OF"
                elif prediction == "CONTAINED_BY":
                    prediction = "CONTAINER_OF"
                elif prediction == "CONTAINER_OF":
                    prediction = "CONTAINED_BY"
                predictions.append((acc_2, acc_1, prediction))

    # Populating METHOD_PREDICTION
    con.drop_table(schema, "METHOD_PREDICTION")
    con.execute(
        """
        CREATE TABLE {}.METHOD_PREDICTION
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            RELATION VARCHAR2(15) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    for i in range(0, len(predictions), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT /*+APPEND*/
            INTO {}.METHOD_PREDICTION (METHOD_AC1, METHOD_AC2, RELATION)
            VALUES (:1, :2, :3)
            """.format(schema),
            predictions[i:i + BULK_INSERT_SIZE]
        )
        con.commit()

    # Optimising METHOD_PREDICTION
    con.execute(
        """
        ALTER TABLE {}.METHOD_PREDICTION
        ADD CONSTRAINT PK_METHOD_PREDICTION
        PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        """.format(schema)
    )
    con.optimise_table(schema, "METHOD_PREDICTION", cascade=True)
    con.grant("SELECT", schema, "METHOD_PREDICTION", "INTERPRO_SELECT")

    # Populating METHOD_MATCH
    con.drop_table(schema, "METHOD_MATCH")
    con.execute(
        """
        CREATE TABLE {}.METHOD_MATCH
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            N_MATCHES NUMBER NOT NULL,
            N_PROT NUMBER NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    signatures = [
        (acc, s["matches"], s["proteins"])
        for acc, s in signatures.items()
    ]
    for i in range(0, len(signatures), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.METHOD_MATCH (
                METHOD_AC, N_MATCHES, N_PROT
            )
            VALUES (:1, :2, :3)
            """.format(schema),
            signatures[i:i + BULK_INSERT_SIZE]
        )
        con.commit()

    # Optimising METHOD_MATCH
    con.execute(
        """
        ALTER TABLE {}.METHOD_MATCH
        ADD CONSTRAINT PK_METHOD_MATCH PRIMARY KEY (METHOD_AC)
        """.format(schema)
    )
    con.optimise_table(schema, "METHOD_MATCH", cascade=True)
    con.grant("SELECT", schema, "METHOD_MATCH", "INTERPRO_SELECT")

    # Creating METHOD_OVERLAP
    con.drop_table(schema, "METHOD_OVERLAP")
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
        """.format(schema)
    )

    # Computing data for METHOD_OVERLAP
    overlaps = []
    for acc_1 in comparisons:
        for acc_2 in comparisons[acc_1]:
            c = comparisons[acc_1][acc_2]

            if c["over"]:
                avg_frac1 = 100 * c["frac_1"] / c["over"]
                avg_frac2 = 100 * c["frac_2"] / c["over"]
                avg_over = c["length"] / c["over"]
            else:
                avg_frac1 = avg_frac2 = avg_over = 0

            """
            Cast to float as cx_Oracle 6.1 throws TypeError
                (expecting integer)
            when the value of the 1st record is an integer
                (e.g. AVG_OVER=0)
            and the value for the second is not (e.g. AVG_OVER=25.6)
            """
            overlaps.append((
                acc_1, acc_2, c["prot"], c["over"], c["prot_over"],
                float(avg_over), float(avg_frac1), float(avg_frac2)
            ))

            if acc_1 != acc_2:
                overlaps.append((
                    acc_2, acc_1, c["prot"], c["over"], c["prot_over"],
                    float(avg_over), float(avg_frac2), float(avg_frac1)
                ))

    # Populating METHOD_OVERLAP
    for i in range(0, len(overlaps), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.METHOD_OVERLAP (
                METHOD_AC1, METHOD_AC2, N_PROT, N_OVER,
                N_PROT_OVER, AVG_OVER, AVG_FRAC1, AVG_FRAC2
            )
            VALUES (:1, :2, :3, :4, :5, :6, :7, :8)
            """.format(schema),
            overlaps[i:i + BULK_INSERT_SIZE]
        )
        con.commit()

    # Optimising METHOD_OVERLAP
    con.execute(
        """
        ALTER TABLE {}.METHOD_OVERLAP
        ADD CONSTRAINT PK_METHOD_OVERLAP
        PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        """.format(schema)
    )
    con.optimise_table(schema, "METHOD_OVERLAP", cascade=True)
    con.grant("SELECT", schema, "METHOD_OVERLAP", "INTERPRO_SELECT")


def calculate_similarities(dsn, schema, coverages, overlaps):
    con.drop_table(schema, "METHOD_SIMILARITY")
    con.execute(
        """
        CREATE TABLE {}.METHOD_SIMILARITY
        (
            METHOD_AC1 VARCHAR2(25) NOT NULL,
            METHOD_AC2 VARCHAR2(25) NOT NULL,
            SIMILARITY BINARY_DOUBLE NOT NULL,
            CONTAINMENT1 BINARY_DOUBLE NOT NULL,
            CONTAINMENT2 BINARY_DOUBLE NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    data = []
    for acc_1 in overlaps:
        for acc_2 in overlaps[acc_1]:
            i = overlaps[acc_1][acc_2]
            u = coverages[acc_1] + coverages[acc_2] - i
            index = i / u
            cont_1 = i / coverages[acc_1]
            cont_2 = i / coverages[acc_2]
            data.append((acc_1, acc_2, index, cont_1, cont_2))

            if len(data) == BULK_INSERT_SIZE:
                con.executemany(
                    """
                    INSERT /*+APPEND*/
                    INTO {}.METHOD_SIMILARITY
                    VALUES (:1, :2, :3, :4, :5)
                    """.format(schema),
                    data
                )
                con.commit()
                data = []

    if data:
        con.executemany(
            """
            INSERT /*+APPEND*/
            INTO {}.METHOD_SIMILARITY
            VALUES (:1, :2, :3, :4, :5)
            """.format(schema),
            data
        )
        con.commit()
        data = []

    con.execute(
        """
        ALTER TABLE {}.METHOD_SIMILARITY
        ADD CONSTRAINT PK_METHOD_SIMILARITY
        PRIMARY KEY (METHOD_AC1, METHOD_AC2)
        """.format(schema)
    )


def optimise_method2protein(dsn, schema):
    con = Connection(dsn)
    con.execute(
        """
        ALTER TABLE {}.METHOD2PROTEIN
        ADD CONSTRAINT PK_METHOD2PROTEIN
        PRIMARY KEY (METHOD_AC, PROTEIN_AC)
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$M
        ON {}.METHOD2PROTEIN (METHOD_AC)
        NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$P
        ON {}.METHOD2PROTEIN (PROTEIN_AC)
        NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$LN
        ON {}.METHOD2PROTEIN (LEFT_NUMBER)
        NOLOGGING
        """.format(schema)
    )
    con.execute(
        """
        CREATE INDEX I_METHOD2PROTEIN$M$DB
        ON {}.METHOD2PROTEIN (METHOD_AC, DBCODE)
        NOLOGGING
        """.format(schema)
    )

    con.optimise_table(schema, "METHOD2PROTEIN", cascade=True)
    con.grant("SELECT", schema, "METHOD2PROTEIN", "INTERPRO_SELECT")


def load_description_counts(dsn, schema, organisers):
    con = Connection(dsn)
    con.drop_table(schema, "METHOD_DESC")
    con.execute(
        """
        CREATE TABLE {}.METHOD_DESC
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            DESC_ID NUMBER(10) NOT NULL,
            REVIEWED_COUNT NUMBER(10) NOT NULL,
            UNREVIEWED_COUNT NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    signatures = {}
    for acc, descriptions in heapq.merge(*organisers, key=lambda x: x[0]):
        if acc in signatures:
            s = signatures[acc]
        else:
            s = signatures[acc] = {}

        for desc_id, dbcode in descriptions:
            if desc_id not in s:
                s[desc_id] = {'S': 0, 'T': 0}

            s[desc_id][dbcode] += 1

    data = []
    for method_acc, descriptions in signatures.items():
        for desc_id, dbcodes in descriptions.items():
            data.append((method_acc, desc_id, dbcodes.get('S', 0),
                         dbcodes.get('T', 0)))

            if len(data) == BULK_INSERT_SIZE:
                con.executemany(
                    """
                    INSERT /*+APPEND*/ INTO {}.METHOD_DESC
                    VALUES (:1, :2, :3, :4)
                    """.format(schema),
                    data
                )
                con.commit()
                data = []

    if data:
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.METHOD_DESC
            VALUES (:1, :2, :3, :4)
            """.format(schema),
            data
        )
        con.commit()

    con.execute(
        """
        ALTER TABLE {}.METHOD_DESC
        ADD CONSTRAINT PK_METHOD_DESC
        PRIMARY KEY (METHOD_AC, DESC_ID)
        NOLOGGING
        """.format(schema)
    )
    con.optimise_table(schema, "METHOD_DESC", cascade=True)
    con.grant("SELECT", schema, "METHOD_DESC", "INTERPRO_SELECT")


def load_taxonomy_counts(dsn, schema, organisers):
    con = Connection(dsn)
    con.drop_table(schema, "METHOD_TAXA")
    con.execute(
        """
        CREATE TABLE {}.METHOD_TAXA
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            RANK VARCHAR2(50) NOT NULL,
            TAX_ID NUMBER(10) NOT NULL,
            PROTEIN_COUNT NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    signatures = {}
    for acc, taxa in heapq.merge(*organisers, key=lambda x: x[0]):
        if acc in signatures:
            s = signatures[acc]
        else:
            s = signatures[acc] = {}

        for rank, tax_id in taxa:
            if rank in s:
                r = s[rank]
            else:
                r = s[rank] = {}

            if tax_id in r:
                r[tax_id] += 1
            else:
                r[tax_id] = 1

    data = []
    for method_acc, ranks in signatures.items():
        for rank, taxa in ranks.items():
            for tax_id, count in taxa.items():
                data.append((method_acc, rank, tax_id, count))

                if len(data) == BULK_INSERT_SIZE:
                    con.executemany(
                        """
                        INSERT /*+APPEND*/ INTO {}.METHOD_TAXA
                        VALUES (:1, :2, :3, :4)
                        """.format(schema),
                        data
                    )
                    con.commit()
                    data = []

    if data:
        con.executemany(
            """
            INSERT /*+APPEND*/ INTO {}.METHOD_TAXA
            VALUES (:1, :2, :3, :4)
            """.format(schema),
            data
        )
        con.commit()

    con.execute(
        """
        ALTER TABLE {}.METHOD_TAXA
        ADD CONSTRAINT PK_METHOD_TAXA
        PRIMARY KEY (METHOD_AC, RANK, TAX_ID)
        NOLOGGING
        """.format(schema)
    )
    con.optimise_table(schema, "METHOD_TAXA", cascade=True)
    con.grant("SELECT", schema, "METHOD_TAXA", "INTERPRO_SELECT")


def process_proteins(dsn, schema, proteins_src, matches_src, processes,
                     **kwargs) -> Tuple[dict, dict, dict, dict, list, list]:
    chunk_size = kwargs.get("chunk_size", 10000)
    dir = kwargs.get("dir")
    max_gap = kwargs.get("max_gap", 20)

    processes = max(1, processes - 1)
    in_queue = Queue(maxsize=processes)
    out_queue = Queue()
    consumers = []
    for _ in range(processes):
        c = ProteinConsumer(dsn, schema, max_gap, in_queue, out_queue,
                            dir=dir)
        c.start()
        consumers.append(c)

    chunk = []
    cnt = 0
    ts = time.time()
    with io.Store(proteins_src) as proteins, io.Store(matches_src) as matches:
        for acc, length, protein_dbcode, desc_id, left_num in proteins:
            _acc, _matches = next(matches)
            while _acc < acc:
                _acc, _matches = next(matches)

            if _acc != acc:
                continue  # Not matches for this proteins

            protein_matches = []
            methods = {}
            methods_fragments = {}

            for m in _matches:
                (method_acc, method_dbcode, method_type,
                 pos_start, pos_end, fragments_str) = m

                if method_dbcode in ('F', 'V'):
                    """
                    PANTHER & PRINTS:
                        Merge protein matches.
                        If the signature is a family*,
                            use the entire protein.

                        * all PANTHER signatures,
                            almost all PRINTS signatures
                    """
                    if method_acc not in methods:
                        methods_fragments[method_acc] = fragments_str

                        if method_type == 'F':
                            # Families: use the entire protein sequence
                            methods[method_acc] = [(1, length)]
                        else:
                            methods[method_acc] = [(pos_start, pos_end)]
                    elif method_type != 'F':
                        # Only for non-families
                        # (because families use the entire protein)
                        methods[method_acc].append((pos_start, pos_end))
                else:
                    protein_matches.append((method_acc, pos_start, pos_end,
                                            fragments_str))

            # Merge matches
            for method_acc in methods:
                min_pos = None
                max_pos = None

                for pos_start, pos_end in methods[method_acc]:
                    if min_pos is None or pos_start < min_pos:
                        min_pos = pos_start
                    if max_pos is None or pos_end > max_pos:
                        max_pos = pos_end

                protein_matches.append((method_acc, min_pos, max_pos,
                                        methods_fragments[method_acc]))

            chunk.append((acc, protein_dbcode, length, desc_id, left_num,
                          protein_matches))

            if len(chunk) == chunk_size:
                in_queue.put(chunk)
                chunk = []

            cnt += 1
            if not cnt % 10000000:
                logging.info("{:>12} ({:.0f} proteins/sec)".format(
                    cnt,
                    cnt / (time.time() - ts)
                ))

    if chunk:
        in_queue.put(chunk)
        chunk = []

    for _ in consumers:
        in_queue.put(None)

    logging.info("{:>12} ({:.0f} proteins/sec)".format(
        cnt,
        cnt / (time.time() - ts)
    ))

    # For predictions
    signatures = {}
    comparisons = {}
    # For residue overlap comparison
    residue_coverages = {}
    residue_overlaps = {}
    # For counts
    name_organisers = []
    taxon_organisers = []
    for _ in consumers:
        s, c, rc, ro, no, to = out_queue.get()

        # Merge dictionaries
        for acc in s:
            if acc in signatures:
                for k, v in s[acc].items():
                    signatures[acc][k] += v
            else:
                signatures[acc] = s[acc]

        for acc_1 in c:
            if acc_1 in comparisons:
                d = comparisons[acc_1]
            else:
                comparisons[acc_1] = c[acc_1]
                continue

            for acc_2 in c[acc_1]:
                if acc_2 in d:
                    for k, v in c[acc_1][acc_2].items():
                        d[acc_2][k] += v
                else:
                    d[acc_2] = c[acc_1][acc_2]

        for acc in rc:
            if acc in residue_coverages:
                residue_coverages[acc] += rc[acc]
            else:
                residue_coverages[acc] = rc[acc]

        for acc_1 in ro:
            if acc_1 in residue_overlaps:
                for acc_2 in ro[acc_1]:
                    if acc_2 in residue_overlaps[acc_1]:
                        residue_overlaps[acc_1][acc_2] += ro[acc_1][acc_2]
                    else:
                        residue_overlaps[acc_1][acc_2] = ro[acc_1][acc_2]
            else:
                residue_overlaps[acc_1] = ro[acc_1]

        name_organisers.append(no)
        taxon_organisers.append(to)

    for c in consumers:
        c.join()

    return (signatures, comparisons, residue_coverages, residue_overlaps,
            name_organisers, taxon_organisers)


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

        for e in sorted(entries.values(), key=_get_entry_key):
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                e["acc"],
                e["type"],
                e["checked"],
                len(e["gained"]),
                len(e["lost"]),
                " | ".join(e["gained"]),
                " | ".join(e["lost"])
            ))


def _get_entry_key(entry):
    if entry["type"] == "F":
        return 0, entry["type"], entry["acc"]
    else:
        return 1, entry["type"], entry["acc"]
