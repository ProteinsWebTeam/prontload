#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.4.3"

import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def cli():
    import argparse
    import os
    import json
    import sys
    from datetime import datetime
    from tempfile import gettempdir

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
    parser.add_argument("-o", "--output", default=default_report,
                        help="output SwissProt report for curators "
                             "(default: {})".format(default_report))
    args = parser.parse_args()

    os.makedirs(args.tmpdir, exist_ok=True)

    with open(args.config, "rt") as fh:
        config = json.load(fh)

    dsn = config["dsn"]
    schema = config["schema"]
    max_gap = int(config["max_gap"])
    
    steps = [
        {
            "name": "clear",
            "func": interpro.clear_schema,
            "args": (dsn, schema),
            "skip": True
        },
        {
            "name": "synonyms",
            "func": interpro.create_synonyms,
            "args": (dsn, "INTERPRO", schema, (
                "ENTRY",
                "ENTRY2METHOD",
                "ENTRY2ENTRY",
                "ENTRY2COMP"
            ))
        },
        {
            "name": "databases",
            "func": interpro.load_databases,
            "args": (dsn, schema)
        },
        {
            "name": "signatures",
            "func": interpro.load_signatures,
            "args": (dsn, schema)
        },
        {
            "name": "taxa",
            "func": interpro.load_taxa,
            "args": (dsn, schema)
        },
        {
            "name": "proteins",
            "func": interpro.load_proteins,
            "args": (dsn, schema)
        },
        {
            "name": "comments",
            "func": uniprot.load_comments,
            "args": (dsn, schema)
        },
        {
            "name": "descriptions",
            "func": uniprot.load_descriptions,
            "args": (dsn, schema, args.tmpdir)
        },
        {
            "name": "enzymes",
            "func": uniprot.load_enzymes,
            "args": (dsn, schema)
        },
        {
            "name": "annotations",
            "func": goa.load_annotations,
            "args": (dsn, schema)
        },
        {
            "name": "publications",
            "func": goa.load_publications,
            "args": (dsn, schema)
        },
        {
            "name": "terms",
            "func": goa.load_terms,
            "args": (dsn, schema)
        },
        {
            # Requires "descriptions" and "taxa"
            "name": "matches",
            "func": interpro.load_matches,
            "args": (dsn, schema),
            "kwargs": dict(
                max_gap=max_gap,
                tmpdir=args.tmpdir
            )
        },
        {
            "name": "copy",
            "func": interpro.copy_schema,
            "args": (dsn, schema)
        },
        {
            "name": "report",
            "func": interpro.report_description_changes,
            "args": (dsn, schema, args.output)
        }
    ]

    step_names = [s["name"] for s in steps]
    to_run = []

    if args.steps:
        for s in args.steps:
            if s in step_names:
                to_run.append(step_names.index(s))
            else:
                sys.stderr.write(
                    "error: invalid step: '{}' "
                    "(choose from {})\n".format(
                        s,
                        ", ".join(["'{}'".format(_s) for _s in step_names])
                    )
                )
                exit(1)
    else:
        to_run = [i for i, s in enumerate(steps) if not s.get("skip")]

    for i in to_run:
        step = steps[i]
        logging.info("running '{}'".format(step["name"]))
        step["func"](*step["args"], **step.get("kwargs", {}))

    logging.info("complete")
