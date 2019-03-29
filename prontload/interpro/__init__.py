import time
from multiprocessing import Process, Queue
from threading import Thread

from .. import get_logger
from ..oracledb import Connection, BULK_INSERT_SIZE

from .utils import enable_schema, get_entry_key

logger = get_logger()


def clear_schema(dsn, schema):
    con = Connection(dsn)
    for table in con.get_tables(schema):
        con.drop_table(schema, table)


def copy_schema(dsn, schema):
    con = Connection(dsn)

    logger.debug("copy                exporting")
    proc = "{}.copy_interpro_analysis.exp_interpro_analysis".format(schema)
    con.exec(proc)
    enable_schema(con, schema)

    logger.debug("copy                importing")
    con.exec("interpro_analysis.drop_all")
    to_drop = []
    for row in con.get(
            """
            SELECT OWNER_NAME, JOB_NAME
            FROM DBA_DATAPUMP_JOBS
            WHERE STATE = 'NOT RUNNING'
            AND ATTACHED_SESSIONS = 0
            """
    ):
        to_drop.append(row)

    for owner, table in to_drop:
        if owner == schema:
            con.drop_table(schema, table)

    proc = "{}.copy_interpro_analysis.imp_interpro_analysis_load".format(schema)
    con.exec(proc)

    enable_schema(con, "INTERPRO_ANALYSIS")


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

        for e in sorted(entries.values(), key=get_entry_key):
            fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                e["acc"],
                e["type"],
                e["checked"],
                len(e["gained"]),
                len(e["lost"]),
                " | ".join(e["gained"]),
                " | ".join(e["lost"])
            ))


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
            IS_READY CHAR(1) DEFAULT 'N',
            CONSTRAINT PK_DATABASE PRIMARY KEY (DBCODE)
        ) NOLOGGING
        """.format(schema)
    )

    con.execute(
        """
        INSERT /*+ APPEND */ INTO {}.CV_DATABASE (
            DBCODE, DBNAME, DBSHORT, VERSION, FILE_DATE
        )
        SELECT DB.DBCODE, DB.DBNAME, DB.DBSHORT, V.VERSION, V.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V
          ON DB.DBCODE = V.DBCODE
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
    logger.debug("matches             inserting protein/feature matches")
    con.execute(
        """
        INSERT /*+ APPEND */ INTO {}.MATCH
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE,
          CASE
            WHEN MODEL_AC IS NOT NULL AND MODEL_AC != METHOD_AC
            THEN MODEL_AC
            ELSE NULL
          END AS MODEL_AC,
          POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.MATCH
        UNION ALL
        SELECT
          PROTEIN_AC, METHOD_AC, DBCODE,
          NULL,
          POS_FROM, POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH
        WHERE DBCODE = 'g'
        """.format(schema)
    )
    con.commit()

    # Finalizing table
    logger.debug("matches             optimising MATCH")
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
    logger.debug("matches             counting #proteins/signature")
    signatures = {}
    for acc, num_proteins in con.get(
            """
            SELECT METHOD_AC, COUNT(DISTINCT PROTEIN_AC)
            FROM {}.MATCH
            GROUP BY METHOD_AC
            """.format(schema)
    ):
        signatures[acc] = num_proteins

    logger.debug("matches             updating METHOD")
    for acc, num_proteins in signatures.items():
        con.execute(
            """
            UPDATE {0}.METHOD
            SET PROTEIN_COUNT = :1
            WHERE METHOD_AC = :2
            """.format(schema),
            num_proteins, acc
        )
    con.commit()


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
        INSERT /*+ APPEND */ INTO {}.PROTEIN (
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
            ABSTRACT VARCHAR2(4000),
            ABSTRACT_LONG CLOB,
            PROTEIN_COUNT NUMBER(8) DEFAULT 0 NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    con.execute(
        """
        INSERT /*+ APPEND */ INTO {}.METHOD (
            METHOD_AC, NAME, DBCODE, CANDIDATE,
            DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
        )
        SELECT
            METHOD_AC, NAME, DBCODE, CANDIDATE,
            DESCRIPTION, SIG_TYPE, ABSTRACT, ABSTRACT_LONG
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
        )
        """.format(schema)
    )

    con.execute(
        """
        INSERT /*+ APPEND */ INTO {}.ETAXI (
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
            INSERT INTO {}.LINEAGE (LEFT_NUMBER, TAX_ID, RANK)
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


def load_method2protein(dsn: str, schema: str, chunk_size: int=10000,
                        dir: str=None, max_gap: int=20, processes: int=1):
    n_workers = max(1, processes - 1)
    task_queue = Queue(maxsize=n_workers)
    done_queue = Queue()
    consumers = []
    for _ in range(n_workers):
        c = utils.ProteinConsumer(dsn, schema, max_gap, task_queue,
                                  done_queue, dir=dir)
        c.start()
        consumers.append(c)

    con = Connection(dsn)
    con.drop_table(schema, "METHOD2PROTEIN")
    con.execute(
        """
        CREATE TABLE {}.METHOD2PROTEIN
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            PROTEIN_AC VARCHAR2(15) NOT NULL,
            DBCODE CHAR(1) NOT NULL,
            MD5 VARCHAR(32) NOT NULL,
            LEN NUMBER(5) NOT NULL,
            LEFT_NUMBER NUMBER NOT NULL,
            DESC_ID NUMBER(10) NOT NULL
        ) NOLOGGING
        """.format(schema)
    )

    query = """
        SELECT
            P.PROTEIN_AC, P.LEN, P.DBCODE, D.DESC_ID,
            NVL(E.LEFT_NUMBER, 0),
            MA.METHOD_AC, MA.POS_FROM, MA.POS_TO, MA.FRAGMENTS
        FROM {0}.PROTEIN P
        INNER JOIN {0}.ETAXI E
          ON P.TAX_ID = E.TAX_ID
        INNER JOIN {0}.PROTEIN_DESC D
          ON P.PROTEIN_AC = D.PROTEIN_AC
        INNER JOIN INTERPRO.MATCH MA
          ON P.PROTEIN_AC = MA.PROTEIN_AC
        WHERE P.FRAGMENT = 'N'
        ORDER BY P.PROTEIN_AC
    """.format(schema)

    _p_acc = None
    p_length = None
    p_dbcode = None
    desc_id = None
    left_num = None
    matches = []
    chunk = []
    cnt = 0
    ts = time.time()
    for row in con.get(query):
        p_acc = row[0]
        if p_acc != _p_acc:
            if _p_acc:
                # Add protein to chunk
                chunk.append((_p_acc, p_dbcode, p_length, desc_id, left_num,
                              matches))

                if len(chunk) == chunk_size:
                    task_queue.put(chunk)
                    chunk = []
            else:
                logger.debug("method2protein      query time: "
                             "{:.0f} seconds".format(time.time()-ts))
                ts = time.time()

            p_length = row[1]
            p_dbcode = row[2]
            desc_id = row[3]
            left_num = row[4]
            matches = []
            _p_acc = p_acc
            cnt += 1
            if not cnt % 10000000:
                logger.debug(
                    "method2protein      "
                    "{:>12} ({:.0f} proteins/sec)".format(
                        cnt, cnt / (time.time() - ts)
                    )
                )

        matches.append((row[5], row[6], row[7], row[8]))

    if _p_acc:
        # Last protein (_p_acc is None only if 0 proteins)
        chunk.append((_p_acc, p_dbcode, p_length, desc_id, left_num,
                      matches))
        task_queue.put(chunk)
        cnt += 1

    for _ in consumers:
        task_queue.put(None)

    logger.debug(
        "method2protein      "
        "{:>12} ({:.0f} proteins/sec)".format(
            cnt, cnt / (time.time() - ts)
        )
    )

    # number of proteins and matches for each signatures
    signatures = {}
    # multiple comparison metrics between two signatures
    comparisons = {}
    # number of residues matched for each signature
    residue_coverages = {}
    # residue overlap between two signatures
    residue_overlaps = {}
    # names organisers
    names = []
    # taxa organisers
    taxa = []
    # temporary space used (in bytes)
    size = 0
    for _ in consumers:
        s, c, rc, ro, n, t, tmp = done_queue.get()
        size += tmp

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

        names.append(n)
        taxa.append(t)

    for c in consumers:
        c.join()

    logger.debug("method2protein      "
                 "temporary disk space: {:,} bytes".format(size))

    # Create the prediction/overlap tables
    utils.make_predictions(con, schema, signatures, comparisons)
    signatures = comparisons = None

    # Calculate similarities
    utils.calculate_similarities(con, schema, residue_coverages, residue_overlaps)
    residue_coverages = residue_overlaps = None

    # Create table (MV-like) for only SwissProt proteins
    utils.create_method2swissprot(dsn, schema)

    t1 = Thread(target=utils.optimise_method2protein, args=(dsn, schema))
    t2 = Thread(target=utils.optimise_method2swissprot, args=(dsn, schema))
    p1 = Process(target=utils.load_description_counts,
                 args=(dsn, schema, names))
    p2 = Process(target=utils.load_taxonomy_counts, args=(dsn, schema, taxa))

    for x in (t1, t2, p1, p2):
        x.start()

    for x in (t1, t2, p1, p2):
        x.join()

    for o in names:
        o.remove()

    for o in taxa:
        o.remove()
