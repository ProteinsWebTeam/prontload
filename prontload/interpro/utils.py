import hashlib
import heapq
from multiprocessing import Process

from .. import get_logger, io
from ..oracledb import Connection, BULK_INSERT_SIZE

logger = get_logger()


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
        signature_info = {}
        accessions = []
        bucket_size = 1000
        cnt = bucket_size
        for acc, dbcode, s_type in con.get(
                """
                SELECT METHOD_AC, DBCODE, SIG_TYPE
                FROM {}.METHOD
                ORDER BY METHOD_AC
                """.format(self.schema)
        ):
            signature_info[acc] = (dbcode, s_type)
            accessions.append(acc)
            if cnt == bucket_size:
                accessions.append(acc)
                cnt = 1
            else:
                cnt += 1

        # Creating organisers
        names = io.Organiser(accessions, dir=self.dir)
        taxa = io.Organiser(accessions, dir=self.dir)

        signatures = {}
        comparisons = {}
        res_coverages = {}
        res_overlaps = {}
        data = []
        for chunk in iter(self.queue_in.get, None):
            for args in chunk:
                protein_acc, dbcode, length, desc_id, left_num, matches = args

                md5 = self.hash(matches, self.max_gap)
                self.compare_fragments(matches, res_coverages, res_overlaps)

                # Merge PANTHER and PRINTS matches, then compare matches
                matches = self.merge_matches(matches, signature_info, length)
                _signatures, _comparisons = self.compare(matches)

                for signature_acc in _signatures:
                    # Taxonomic origins
                    taxa.add(signature_acc, left_num)

                    # UniProt descriptions
                    names.add(signature_acc, (desc_id, dbcode))

                    data.append((
                        signature_acc, protein_acc, dbcode,
                        md5, length, left_num, desc_id
                    ))

                    if signature_acc not in signatures:
                        signatures[signature_acc] = {
                            "proteins": 0,
                            "matches": 0
                        }

                    signatures[signature_acc]["proteins"] += 1
                    n_matches = _signatures[signature_acc]
                    signatures[signature_acc]["matches"] += n_matches

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

            names.dump()
            taxa.dump()

            while len(data) >= BULK_INSERT_SIZE:
                con.executemany(
                    """
                    INSERT INTO {}.METHOD2PROTEIN (
                        METHOD_AC, PROTEIN_AC, DBCODE,
                        MD5, LEN, LEFT_NUMBER, DESC_ID
                    )
                    VALUES (:1, :2, :3, :4, :5, :6, :7)
                    """.format(self.schema),
                    data[:BULK_INSERT_SIZE]
                )
                con.commit()
                data = data[BULK_INSERT_SIZE:]

        for i in range(0, len(data), BULK_INSERT_SIZE):
            con.executemany(
                """
                INSERT INTO {}.METHOD2PROTEIN (
                    METHOD_AC, PROTEIN_AC, DBCODE,
                    MD5, LEN, LEFT_NUMBER, DESC_ID
                )
                VALUES (:1, :2, :3, :4, :5, :6, :7)
                """.format(self.schema),
                data[i:i+BULK_INSERT_SIZE]
            )
            con.commit()

        size = names.merge() + taxa.merge()
        self.queue_out.put((
            signatures,
            comparisons,
            res_coverages,
            res_overlaps,
            names,
            taxa,
            size
        ))

    @staticmethod
    def merge_matches(matches: list, signatures: dict, length: int) -> list:
        result = []
        to_merge = {}

        for acc, pos_start, pos_end, _ in matches:
            dbcode, _type = signatures[acc]

            if dbcode in ('F', 'V'):
                """
                PANTHER & PRINTS: Merge protein matches.
                If the signature is a family*, use the entire protein.
                * all PANTHER signatures, almost all PRINTS signatures
                """
                if acc not in to_merge:
                    if _type == 'F':
                        # Families: use the entire protein sequence
                        to_merge[acc] = [(1, length)]
                    else:
                        to_merge[acc] = [(pos_start, pos_end)]
                elif _type != 'F':
                    # Only for non-families
                    # (because families use the entire protein)
                    to_merge[acc].append((pos_start, pos_end))
            else:
                result.append((acc, pos_start, pos_end))

        for acc in to_merge:
            min_pos = None
            max_pos = None

            for pos_start, pos_end in to_merge[acc]:
                if min_pos is None or pos_start < min_pos:
                    min_pos = pos_start

                if max_pos is None or pos_end > max_pos:
                    max_pos = pos_end

            result.append((acc, min_pos, max_pos))

        return result

    @staticmethod
    def compare(matches: list) -> tuple:
        comparisons = {}
        signatures = {}

        for m1_acc, m1_start, m1_end in matches:
            m1_len = m1_end - m1_start + 1

            if m1_acc in signatures:
                signatures[m1_acc] += 1
                m1 = comparisons[m1_acc]
            else:
                signatures[m1_acc] = 1
                m1 = comparisons[m1_acc] = {}

            for m2_acc, m2_start, m2_end in matches:
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
    def hash(matches: list, max_gap: int) -> str:
        # flatten matches
        locations = []

        # Fourth item is fragments_str, which is not used here
        for acc, pos_start, pos_end, _ in matches:
            locations.append((pos_start, acc))
            locations.append((pos_end, acc))

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
        signatures = []

        for pos, acc in locations:
            if pos > offset + max_gap:
                for _pos, _ac in signatures:
                    structure.append(_ac)

                signatures = []
                structure.append('')  # add a gap

            offset = pos
            signatures.append((pos, acc))

        for _pos, _ac in signatures:
            structure.append(_ac)

        return hashlib.md5(
            '/'.join(structure).encode("utf-8")
        ).hexdigest()

    @staticmethod
    def compare_fragments(matches: list, residues: dict, overlaps: dict):
        # Parse the fragments string
        signatures = {}
        for acc, start, end, fragments_str in matches:
            if acc in signatures:
                fragments = signatures[acc]
            else:
                fragments = signatures[acc] = []

            if fragments_str is None:
                fragments.append((start, end))
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
                        fragments.append((pos_start, pos_end))

                if fallback:
                    # Fallback to match start/end positions
                    fragments.append((start, end))

        # Sort fragments by position
        for fragments in signatures.values():
            fragments.sort()

        # Evaluate overlap between all signatures
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


def make_predictions(con, schema, signatures, comparisons):
    logger.debug("method2protein      making predictions")
    candidates = set()
    non_prosite_candidates = set()
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
        )
        """.format(schema)
    )

    for i in range(0, len(predictions), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT
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
        )
        """.format(schema)
    )

    signatures = [
        (acc, s["matches"], s["proteins"])
        for acc, s in signatures.items()
    ]
    for i in range(0, len(signatures), BULK_INSERT_SIZE):
        con.executemany(
            """
            INSERT INTO {}.METHOD_MATCH (
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
        )
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
            INSERT INTO {}.METHOD_OVERLAP (
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
    logger.debug("method2protein      predictions done")


def calculate_similarities(con, schema, coverages, overlaps):
    logger.debug("method2protein      calculating similarities")
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
        )
        """.format(schema)
    )

    data = []
    for acc_1 in overlaps:
        for acc_2 in overlaps[acc_1]:
            i = overlaps[acc_1][acc_2]
            if not i:
                # Not a single residue of overlap
                continue

            u = coverages[acc_1] + coverages[acc_2] - i
            index = i / u
            cont_1 = i / coverages[acc_1]
            cont_2 = i / coverages[acc_2]
            data.append((acc_1, acc_2, index, cont_1, cont_2))

            if len(data) == BULK_INSERT_SIZE:
                con.executemany(
                    """
                    INSERT
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
            INSERT
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
    con.optimise_table(schema, "METHOD_SIMILARITY", cascade=True)
    con.grant("SELECT", schema, "METHOD_SIMILARITY", "INTERPRO_SELECT")
    logger.debug("method2protein      similarities done")


def optimise_method2protein(dsn, schema):
    logger.debug("method2protein      optimising METHOD2PROTEIN")
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

    con.optimise_table(schema, "METHOD2PROTEIN", cascade=True)
    con.grant("SELECT", schema, "METHOD2PROTEIN", "INTERPRO_SELECT")
    logger.debug("method2protein      METHOD2PROTEIN ready")


def create_method2swissprot(dsn, schema):
    logger.debug("method2protein      creating METHOD2SWISSPROT")
    con = Connection(dsn)
    con.drop_table(schema, "METHOD2SWISSPROT")
    con.execute(
        """
        CREATE TABLE {0}.METHOD2SWISSPROT
        AS
        SELECT
          METHOD_AC, PROTEIN_AC, MD5, LEN, LEFT_NUMBER, DESC_ID
        FROM {0}.METHOD2PROTEIN
        WHERE DBCODE = 'S'
        """.format(schema)
    )
    logger.debug("method2protein      METHOD2SWISSPROT created")


def optimise_method2swissprot(dsn, schema):
    logger.debug("method2protein      optimising METHOD2SWISSPROT")
    con = Connection(dsn)
    con.execute(
        """
        ALTER TABLE {}.METHOD2SWISSPROT
        ADD CONSTRAINT PK_METHOD2SWISSPROT
        PRIMARY KEY (METHOD_AC, PROTEIN_AC)
        """.format(schema)
    )
    con.optimise_table(schema, "METHOD2SWISSPROT", cascade=True)
    con.grant("SELECT", schema, "METHOD2SWISSPROT", "INTERPRO_SELECT")
    logger.debug("method2protein      METHOD2SWISSPROT ready")


def enable_schema(con, schema):
    con.execute(
        """
        UPDATE {}.CV_DATABASE
        SET IS_READY = 'Y'
        """.format(schema)
    )
    con.commit()


def get_entry_key(entry):
    if entry["type"] == "F":
        return 0, entry["type"], entry["acc"]
    else:
        return 1, entry["type"], entry["acc"]


def load_description_counts(dsn: str, schema: str, organisers: list):
    logger.debug("method2protein      creating METHOD_DESC")

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
        )
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
    for acc, descriptions in signatures.items():
        for desc_id, dbcodes in descriptions.items():
            data.append((acc, desc_id, dbcodes['S'], dbcodes['T']))

            if len(data) == BULK_INSERT_SIZE:
                con.executemany(
                    """
                    INSERT INTO {}.METHOD_DESC
                    VALUES (:1, :2, :3, :4)
                    """.format(schema),
                    data
                )
                con.commit()
                data = []

    if data:
        con.executemany(
            """
            INSERT INTO {}.METHOD_DESC
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
    logger.debug("method2protein      METHOD_DESC ready")


def load_taxonomy_counts(dsn: str, schema: str, organisers: list):
    logger.debug("method2protein      creating METHOD_TAXA")

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
        )
        """.format(schema)
    )

    # Get lineages for the METHOD_TAXA table
    lineages = {}
    for left_num, tax_id, rank in con.get(
            """
            SELECT LEFT_NUMBER, TAX_ID, RANK
            FROM {}.LINEAGE
            WHERE RANK IN (
              'superkingdom', 'kingdom', 'phylum', 'class', 'order',
              'family', 'genus', 'species'
            )
            """.format(schema)
    ):
        if left_num in lineages:
            lineages[left_num][rank] = tax_id
        else:
            lineages[left_num] = {rank: tax_id}

    signatures = {}
    _ranks = ["superkingdom", "kingdom", "phylum", "class", "order",
              "family", "genus", "species"]
    for acc, left_numbers in heapq.merge(*organisers, key=lambda x: x[0]):
        if acc in signatures:
            s = signatures[acc]
        else:
            s = signatures[acc] = {}

        for left_num in left_numbers:
            ranks = lineages.get(left_num, {"no rank": -1})

            for rank in _ranks:
                tax_id = ranks.get(rank, -1)
                if rank in s:
                    r = s[rank]
                else:
                    r = s[rank] = {}

                if tax_id in r:
                    r[tax_id] += 1
                else:
                    r[tax_id] = 1

    lineages = None
    data = []
    for method_acc, ranks in signatures.items():
        for rank, taxa in ranks.items():
            for tax_id, count in taxa.items():
                data.append((method_acc, rank, tax_id, count))

                if len(data) == BULK_INSERT_SIZE:
                    con.executemany(
                        """
                        INSERT INTO {}.METHOD_TAXA
                        VALUES (:1, :2, :3, :4)
                        """.format(schema),
                        data
                    )
                    con.commit()
                    data = []

    if data:
        con.executemany(
            """
            INSERT INTO {}.METHOD_TAXA
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
    logger.debug("method2protein      METHOD_TAXA ready")
