#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.4.5"

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def exec_functions(*iterable):
    for name, func, args, kwargs in iterable:
        logging.info("{:<20}running".format(name))
        func(*args, **kwargs)
        logging.info("{:<20}done".format(name))


def exec_function(queue):
    while True:
        args = queue.get()
        if args is None:
            break

        func, args = args
        func(*args)


def cli():
    import argparse
    import multiprocessing as mp
    import os
    import json
    from datetime import datetime
    from tempfile import gettempdir
    from threading import Thread

    from . import goa, interpro, uniprot

    default_report = "swissprot_report_{}.tsv".format(
        datetime.today().strftime("%Y_%m_%d")
    )

    parser = argparse.ArgumentParser(
        description="Refresh Pronto with the latest "
                    "data from InterPro, GOA, and UniProt"
    )
    parser.add_argument("config",
                        help="config JSON file")
    parser.add_argument("-s", "--steps", metavar="step", nargs="+",
                        help="steps to perform (default: all)")
    parser.add_argument("-t", "--tmpdir", default=gettempdir(),
                        help="temporary directory "
                             "(default: {})".format(gettempdir()))
    parser.add_argument("-p", "--processes", type=int, default=1,
                        help="number of processes (default: 1)")
    parser.add_argument("-o", "--output", default=default_report,
                        help="output SwissProt report for curators "
                             "(default: {})".format(default_report))
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__),
                        help="show the version and quit")
    args = parser.parse_args()

    os.makedirs(args.tmpdir, exist_ok=True)

    with open(args.config, "rt") as fh:
        config = json.load(fh)

    dsn = config["dsn"]
    schema = config["schema"]
    max_gap = int(config["max_gap"])

    steps = {
        # In priority, if called
        "clear": {
            "func": interpro.clear_schema,
            "args": (dsn, schema),
            "skip": True
        },

        # Grouped in a separate thread
        "databases": {
            "func": interpro.load_databases,
            "args": (dsn, schema),
        },
        "comments": {
            "func": uniprot.load_comments,
            "args": (dsn, schema),
        },
        "enzymes": {
            "func": uniprot.load_enzymes,
            "args": (dsn, schema),
        },
        "annotations": {
            "func": goa.load_annotations,
            "args": (dsn, schema),
        },
        "publications": {
            "func": goa.load_publications,
            "args": (dsn, schema),
        },
        "terms": {
            "func": goa.load_terms,
            "args": (dsn, schema),
        },
        "matches": {
            "func": interpro.load_matches,
            "args": (dsn, schema),
        },

        # In main thread, one after the other
        "synonyms": {
            "func": interpro.create_synonyms,
            "args": (dsn, "INTERPRO", schema)
        },
        "signatures": {
            "func": interpro.load_signatures,
            "args": (dsn, schema),
        },
        "taxa": {
            "func": interpro.load_taxa,
            "args": (dsn, schema),
        },

        # In a separate thread
        "proteins": {
            "func": interpro.load_proteins,
            "args": (dsn, schema),
        },

        # In a separate process
        "descriptions": {
            "func": uniprot.load_descriptions,
            "args": (dsn, schema, args.tmpdir),
        },

        # In the main thread
        "predictions": {
            "func": interpro.load_matches,
            "args": (dsn, schema),
            "kwargs": dict(
                processes=args.processes,
                max_gap=max_gap,
                tmpdir=args.tmpdir
            ),
        },

        # In the main thread
        "report": {
            "func": interpro.report_description_changes,
            "args": (dsn, schema, args.output)
        },

        # After everything is done
        "copy": {
            "func": interpro.copy_schema,
            "args": (dsn, schema)
        }
    }

    if args.steps:
        for step in steps.values():
            step["run"] = False

        for name in args.steps:
            try:
                step = steps[name]
            except KeyError:
                parser.error(
                    "invalid step: '{}' "
                    "(choose from {})\n".format(
                        name,
                        ", ".join(["'{}'".format(s) for s in sorted(steps)])
                    )
                )
            else:
                step["run"] = True
    else:
        for step in steps.values():
            step["run"] = not step.get("skip", False)

    if steps["clear"]["run"]:
        logging.info("{:<20}running".format("clear"))
        s = steps["clear"]
        s["func"](*s["args"], **s.get("kwargs", {}))
        logging.info("{:<20}done".format("clear"))

    # Start background steps in separate thread
    group = []
    for name in ("databases", "comments", "enzymes", "annotations",
                 "publications", "terms", "matches"):
        s = steps[name]
        if s["run"]:
            group.append((name, s["func"], s["args"], s.get("kwargs", {})))

    t_group = Thread(target=exec_functions, args=group)
    t_group.start()

    for name in ("synonyms", "signatures", "taxa"):
        s = steps[name]
        if s["run"]:
            logging.info("{:<20}running".format(name))
            s["func"](*s["args"], **s.get("kwargs", {}))
            logging.info("{:<20}done".format(name))

    s = steps["proteins"]
    if s["run"]:
        logging.info("{:<20}running".format("proteins"))
        t_proteins = Thread(target=s["func"], args=s["args"],
                            kwargs=s.get("kwargs", {}))
        t_proteins.start()
    else:
        t_proteins = None

    s = steps["descriptions"]
    if s["run"]:
        logging.info("{:<20}running".format("descriptions"))
        p_descriptions = mp.Process(target=s["func"], args=s["args"],
                                    kwargs=s.get("kwargs", {}))
        p_descriptions.start()
    else:
        p_descriptions = None

    # Wait until proteins are loaded in Oracle
    if t_proteins:
        t_proteins.join()
        logging.info("{:<20}done".format("proteins"))

    s = steps["predictions"]
    if s["run"]:
        matches_f = os.path.join(args.tmpdir, "matches")
        proteins_f = os.path.join(args.tmpdir, "proteins")

        """
        Export matches in a separate thread
        One process is reserved to descriptions load and proteins export
        """
        logging.info("exporting matches")
        t_matches = Thread(target=interpro.dump_matches,
                           args=(dsn, schema, args.processes-1, matches_f,
                                 args.tmpdir))
        t_matches.start()

        if p_descriptions:
            # Then wait until descriptions are loaded
            p_descriptions.join()
            logging.info("{:<20}done".format("descriptions"))

        # When descriptions are loaded, we can export proteins
        logging.info("exporting proteins")
        p_proteins = mp.Process(target=interpro.dump_proteins,
                                args=(dsn, schema, proteins_f))
        p_proteins.start()

        p_proteins.join()
        t_matches.join()

        # Once both exports are completed, we can process proteins
        logging.info("processing proteins")
        res = interpro.process_proteins(dsn, schema, proteins_f, matches_f,
                                        args.processes, dir=args.tmpdir,
                                        max_gap=max_gap)
        """
        s: dict, number of proteins and matches for each signatures
        c: dict (of dict), multiple comparison metrics between two signatures
        rc: dict, number of residues matched for each signature
        ro: dict (dict), residue overlap between two signatures
        no: list, organisers of UniProt descriptions per signature
        to: list, organisers of taxa/ranks per signature
        """
        s, c, rc, ro, no, to = res

        # We don't need the stores any more
        os.remove(proteins_f)
        os.remove(matches_f)

        # Finalise the METHOD2PROTEIN table in a separate thread
        logging.info("optimising METHOD2PROTEIN")
        t_method2proteins = Thread(target=interpro.optimise_method2protein,
                                   args=(dsn, schema))
        t_method2proteins.start()

        """
        We still have four steps to do:
            - make and store predictions (`s` and `c`)
            - calculate the Jaccard/containment indices (`rc`, `ro`)
            - merge UniProt description organisers and store counts (`no`)
            - merge taxa/ranks organisers and store counts (`to`)

        -> create a pool of workers and submit functions and arguments
        """

        # Pool of between 1 and 4 workers (inclusive)
        n = min(4, max(1, args.processes-1))
        pool = []
        q = mp.Queue(maxsize=n)
        for _ in range(n):
            p = mp.Process(target=exec_function, args=(q,))
            p.start()
            pool.append(p)

        logging.info("making predictions")
        q.put((interpro.make_predictions, (dsn, schema, s, c)))

        """
        On the following calls:
            log message only after the task is accepted in the queue
            so we never log something way before it actually starts
        """
        q.put((interpro.calculate_similarities, (dsn, schema, rc, ro)))
        logging.info("calculating similarities")

        q.put((interpro.load_description_counts, (dsn, schema, no)))
        logging.info("loading UniProt description counts")

        q.put((interpro.load_taxonomy_counts, (dsn, schema, ro)))
        logging.info("loading taxonomic origin counts")

        for _ in pool:
            q.put(None)

        for p in pool:
            p.join()

        t_method2proteins.join()
    elif p_descriptions:
        p_descriptions.join()
        logging.info("{:<20}done".format("descriptions"))

    t_group.join()

    for name in ("report", "copy"):
        s = steps[name]
        if s["run"]:
            logging.info("{:<20}running".format(name))
            s["func"](*s["args"], **s.get("kwargs", {}))
            logging.info("{:<20}done".format(name))

    logging.info("complete")
