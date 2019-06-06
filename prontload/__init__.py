#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.5.6"

import logging


def get_logger(name: str="prontload",
               level: int=logging.INFO) -> logging.Logger:

    logger = logging.getLogger(name)
    if not logger.handlers:
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(fmt="%(asctime)s: %(message)s",
                                          datefmt="%Y-%m-%d %H:%M:%S"))
        logger.addHandler(ch)
        logger.setLevel(level)
    elif logger.level != level:
        logger.setLevel(level)

    return logger


def exec_functions(*args):
    logger = get_logger()

    for name, func, f_args, f_kwargs in args:
        logger.info("{:<20}running".format(name))
        func(*f_args, **f_kwargs)
        logger.info("{:<20}done".format(name))


def cli():
    import argparse
    import os
    import json
    from datetime import datetime
    from multiprocessing import Process
    from tempfile import gettempdir
    from threading import Thread

    from . import goa, interpro, oracledb, uniprot

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
    parser.add_argument("--verbose",
                        help="display additional logging messages",
                        action="store_const", const=logging.DEBUG,
                        default=logging.INFO)
    args = parser.parse_args()

    os.makedirs(args.tmpdir, exist_ok=True)
    get_logger(level=args.verbose)

    with open(args.config, "rt") as fh:
        config = json.load(fh)

    users = config["users"]
    dsn = users["main"] + '@' + config["dsn"]
    schema = users["main"].split('/')[0]

    dsn_dst = users["secondary"] + '@' + config["dsn"]
    schema_dst = users["secondary"].split('/')[0]

    max_gap = int(config["max_gap"])

    steps = {
        # In priority, if called
        "clear": {
            "func": oracledb.clear_schema,
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
        "method2protein": {
            "func": interpro.load_method2protein,
            "args": (dsn, schema),
            "kwargs": {
                "dir": args.tmpdir,
                "max_gap": max_gap,
                "processes": args.processes
            }
        },

        # In the main thread
        "report": {
            "func": interpro.report_description_changes,
            "args": (dsn, schema, args.output)
        },

        "copy": {
            "func": interpro.copy_schema,
            "args": (dsn, schema, dsn_dst, schema_dst)
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

    s = steps["clear"]
    if s["run"]:
        exec_functions(("clear", s["func"], s["args"], s.get("kwargs", {})))

    # Start background steps in separate thread
    group = []
    for name in ("databases", "comments", "enzymes", "annotations",
                 "publications", "terms", "matches"):
        s = steps[name]
        if s["run"]:
            group.append((name, s["func"], s["args"], s.get("kwargs", {})))

    t_group = Thread(target=exec_functions, args=group)
    t_group.start()

    group = []
    for name in ("signatures", "taxa"):
        s = steps[name]
        if s["run"]:
            group.append((name, s["func"], s["args"], s.get("kwargs", {})))
    exec_functions(*group)

    s = steps["proteins"]
    if s["run"]:
        _args = ("proteins", s["func"], s["args"], s.get("kwargs", {}))
        t_proteins = Thread(target=exec_functions, args=(_args,))
        t_proteins.start()
    else:
        t_proteins = None

    s = steps["descriptions"]
    if s["run"]:
        _args = ("descriptions", s["func"], s["args"], s.get("kwargs", {}))
        p_descriptions = Process(target=exec_functions, args=(_args,))
        p_descriptions.start()
    else:
        p_descriptions = None

    # Wait until proteins are loaded in Oracle
    if t_proteins:
        t_proteins.join()

    s = steps["method2protein"]
    if s["run"]:
        if p_descriptions:
            # Wait until descriptions are loaded
            p_descriptions.join()

        exec_functions(("method2protein", s["func"], s["args"],
                        s.get("kwargs", {})))
    elif p_descriptions:
        p_descriptions.join()

    t_group.join()

    group = []
    for name in ("report", "copy"):
        s = steps[name]
        if s["run"]:
            group.append((name, s["func"], s["args"], s.get("kwargs", {})))
    exec_functions(*group)

    logging.info("complete")
