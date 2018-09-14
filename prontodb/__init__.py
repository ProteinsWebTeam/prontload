#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.4.0"


def cli():
    import argparse
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
    parser.add_argument("-p", "--threads", default=1, type=int,
                        help="number of threads")
    parser.add_argument("-t", "--temp", default=gettempdir(),
                        help="temporary directory "
                             "(default: {})".format(gettempdir()))
    parser.add_argument("-o", "--output", default=default_report,
                        help="output SwissProt report for curators "
                             "(default: {})".format(default_report))
    args = parser.parse_args()

    with open(args.config, 'rt') as fh:
        config = json.load(fh)

    dsn = config['dsn']
    schema = config['schema']
    max_gap = int(config['max_gap'])
    
    steps = [
        {
            "name": "clear",
            "func": interpro.clear_schema,
            "args": (dsn, schema),
            "skip": True
        },
        {
            "name": "synomyms",
            "func": interpro.create_synonyms,
            "args": (dsn, 'INTERPRO', schema, (
                "ENTRY",
                "ENTRY2METHOD",
                "ENTRY2ENTRY",
                "ENTRY2COMP",
                "METHOD2SWISS_DE"
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
            "args": (dsn, schema, args.temp)
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
            "name": "matches",
            "func": interpro.load_matches,
            "args": (dsn, schema, args.threads, max_gap, args.temp)
        },
        {
            "name": "copy",
            "func": interpro.copy_schema,
            "args": (dsn, schema)
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
                        ', '.join(["'{}'".format(_s) for _s in step_names])
                    )
                )
                exit(1)
    else:
        to_run = [i for i, s in enumerate(steps) if not s.get("skip")]

    for i in to_run:
        step = steps[i]
        sys.stderr.write(
            "{:%y-%m-%d %H:%M:%S}: running '{}'\n".format(
                datetime.now(), step["name"]
            )
        )
        step["func"](*step["args"])

    sys.stderr.write(
        "{:%y-%m-%d %H:%M:%S}: complete\n".format(datetime.now()))

