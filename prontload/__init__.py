#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .interpro import database


_STEPS = [
    # Format: (name, function)

    ('databases', database.load_databases),
    ('synonyms', database.create_synonyms),
    ('comments', database.load_comments),
    ('descriptions', database.load_descriptions),
    ('enzymes', database.load_enzymes),
    ('protein2go', database.load_protein2go),
    ('methods', database.load_methods),
    ('taxonomies', database.load_taxonomies),
    ('terms', database.load_terms),
    ('proteins', database.load_proteins),
    ('matches', database.process_matches),
    ('report', database.report_swiss_descriptions)
]


def get_steps():
    return [step for step, fn in _STEPS]


def run(steps, dsn, schema, **kwargs):
    all_steps = get_steps()

    if steps:
        # Reorder steps
        steps.sort(key=lambda x: all_steps.index(x))
    else:
        steps = all_steps

    all_steps = dict(_STEPS)
    for step in steps:
        fn = all_steps[step]
        fn(dsn, schema, **kwargs)
