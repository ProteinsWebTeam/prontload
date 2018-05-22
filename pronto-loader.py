#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
from tempfile import gettempdir
from datetime import datetime

import prontload


def main():
    default_report = 'swissprot_report_{}.tsv'.format(datetime.today().strftime('%Y_%m_%d'))

    parser = argparse.ArgumentParser(
        description='Refresh Pronto with the latest data from InterPro, GOA, and UniProt'
    )
    parser.add_argument('config',
                        help='config JSON file')
    parser.add_argument('-s', '--steps', metavar='step', nargs='+', choices=prontload.get_steps(),
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

    prontload.run(args.steps, dsn, schema, **kwargs)


if __name__ == '__main__':
    main()
