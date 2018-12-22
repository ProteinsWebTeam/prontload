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


def cli():
    import argparse
    import os
    import json
    from datetime import datetime
    from multiprocessing import Process
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
        "predictions": {},

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
        p_descriptions = Process(target=s["func"], args=s["args"],
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
        if p_descriptions:
            # Wait until descriptions are loaded
            p_descriptions.join()
            logging.info("{:<20}done".format("descriptions"))

        # Create METHOD2PROTEIN table
        logging.info("{:<20}running".format("predictions"))
        res = interpro.load_method2protein(dsn, schema,
                                           dir=args.tmpdir,
                                           max_gap=max_gap,
                                           processes=args.processes)
        logging.info("{:<20}done".format("predictions"))

        """
        s: dict, number of proteins and matches for each signatures
        c: dict (of dict), multiple comparison metrics between two signatures
        rc: dict, number of residues matched for each signature
        ro: dict (dict), residue overlap between two signatures
        """
        s, c, rc, ro = res

        # Finalise the METHOD2PROTEIN table in a separate thread
        logging.info("optimising METHOD2PROTEIN")
        t_method2proteins = Thread(target=interpro.optimise_method2protein,
                                   args=(dsn, schema))
        t_method2proteins.start()

        logging.info("making predictions")
        p1 = Process(target=interpro.make_predictions,
                     args=(dsn, schema, s, c))
        p1.start()

        logging.info("calculating similarities")
        p2 = Process(target=interpro.calculate_similarities,
                     args=(dsn, schema, rc, ro))
        p2.start()

        logging.info("loading count tables")
        interpro.load_count_tables(dsn, schema,
                                   dir=args.tmpdir,
                                   processes=args.processes-2)

        t_method2proteins.join()
        p1.join()
        p2.join()
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
