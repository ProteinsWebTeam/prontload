#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = '0.3'


def cli():
    import argparse
    import json
    import logging
    from datetime import datetime
    from tempfile import gettempdir

    from prontodb import annotations, extra, proteins, signatures

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s: %(levelname)s: %(message)s',
        datefmt='%y-%m-%d %H:%M:%S'
    )

    steps = [
        # Format: (name, function)

        ('databases', extra.load_databases),
        ('synonyms', extra.create_synonyms),
        ('comments', proteins.load_comments),
        ('descriptions', proteins.load_descriptions),
        ('enzymes', proteins.load_enzymes),
        ('protein2go', annotations.load_protein2go),
        ('methods', signatures.load_methods),
        ('taxonomies', extra.load_taxonomies),
        ('terms', annotations.load_terms),
        ('proteins', proteins.load_proteins),
        ('matches', signatures.load_matches_and_predictions),
        ('report', signatures.report_swiss_descriptions)
    ]
    step_names = [step for step, fn in steps]

    default_report = 'swissprot_report_{}.tsv'.format(datetime.today().strftime('%Y_%m_%d'))

    parser = argparse.ArgumentParser(
        description='Refresh Pronto with the latest data from InterPro, GOA, and UniProt'
    )
    parser.add_argument('config',
                        help='config JSON file')
    parser.add_argument('-s', '--steps', metavar='step', nargs='+', choices=step_names,
                        help='steps to perform (default: all)')
    parser.add_argument('-p', '--threads', default=1, type=int,
                        help='number of threads (\'matches\' step only)')
    parser.add_argument('-t', '--temp', default=gettempdir(),
                        help='temporary directory (default: {})'.format(gettempdir()))
    parser.add_argument('-o', '--output', default=default_report,
                        help='output SwissProt report for curators (default: {})'.format(default_report))
    args = parser.parse_args()

    with open(args.config, 'rt') as fh:
        config = json.load(fh)

    dsn = config['dsn']
    schema = config['schema']
    max_gap = int(config['max_gap'])

    kwargs = {
        'chunk_size': 100000,
        'tmpdir': args.temp,
        'max_gap': max_gap,
        'processes': args.threads,
        'output': args.output
    }

    to_run = args.steps

    if to_run:
        # Reorder steps
        to_run.sort(key=lambda x: step_names.index(x))
    else:
        to_run = step_names

    steps = dict(steps)
    for step in to_run:
        fn = steps[step]
        fn(dsn, schema, **kwargs)
