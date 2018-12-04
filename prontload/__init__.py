#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.4.5"

import logging
from typing import Iterable

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def exec_functions(*iterable: Iterable):
    for name, func, args, kwargs in iterable:
        logging.info("running '{}'".format(name))
        func(*args, **kwargs)


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
        logging.info("running 'clear'")
        s = steps["clear"]
        s["func"](*s["args"], **s.get("kwargs", {}))

    # Start background steps in separate thread
    group = []
    for name in ("databases", "comments", "enzymes", "annotations",
                 "publications", "terms", "matches"):
        s = steps[name]
        group.append((name, s["func"], s["args"], s.get("kwargs", {})))

    t_group = Thread(target=exec_functions, args=group)
    t_group.start()

    for name in ("synonyms", "signatures", "taxa"):
        s = steps[name]
        if s["run"]:
            logging.info("running '{}'".format(name))
            s["func"](*s["args"], **s.get("kwargs", {}))

    s = steps["proteins"]
    if s["run"]:
        logging.info("running 'proteins'")
        t_proteins = Thread(target=s["func"], args=s["args"],
                            kwargs=s.get("kwargs", {}))
    else:
        t_proteins = Thread()
    t_proteins.start()

    s = steps["descriptions"]
    if s["run"]:
        logging.info("running 'descriptions'")
        p_descriptions = Process(target=s["func"], args=s["args"],
                                 kwargs=s.get("kwargs", {}))
    else:
        p_descriptions = Process()
    p_descriptions.start()

    # Wait until proteins are loaded in Oracle
    t_proteins.join()

    s = steps["predictions"]
    if s["run"]:
        matches_f = os.path.join(args.tmpdir, "matches")
        proteins_f = os.path.join(args.tmpdir, "proteins")

        """
        Export matches in a separate thread
        One process is reserved to descriptions load and proteins export
        """
        logging.info("exporting matches")
        processes = args.processes - 1
        t_matches = Thread(target=interpro.dump_matches,
                           args=(dsn, schema, processes, matches_f,
                                 args.tmpdir))
        t_matches.start()

        # Then wait until descriptions are loaded
        p_descriptions.join()

        # When descriptions are loaded, we can export proteins
        logging.info("exporting proteins")
        p_proteins = Process(target=interpro.dump_proteins,
                             args=(dsn, schema, proteins_f))
        p_proteins.start()

        p_proteins.join()
        t_matches.join()
    else:
        p_descriptions.join()

    t_group.join()

    logging.info("complete")
