#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os


from . import oracledb


def load_methods(dsn, schema, **kwargs):
    logging.info('loading member database signatures')

    con = oracledb.connect(dsn)
    cur = con.cursor()
    oracledb.drop_table(cur, schema, 'METHOD')

    cur.execute(
        """
        CREATE TABLE {}.METHOD
        (
            METHOD_AC VARCHAR2(25) NOT NULL,
            NAME VARCHAR2(100),
            DBCODE CHAR(1) NOT NULL,
            CANDIDATE CHAR(1) NOT NULL,
            DESCRIPTION VARCHAR2(220),
            SIG_TYPE CHAR(1)
        ) NOLOGGING
        """.format(schema)
    )

    cur.execute(
        """
        INSERT /*+APPEND*/ INTO {}.METHOD (METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE)
        SELECT METHOD_AC, NAME, DBCODE, CANDIDATE, DESCRIPTION, SIG_TYPE
        FROM INTERPRO.METHOD
        """.format(schema)
    )
    con.commit()

    cur.execute(
        """
        ALTER TABLE {}.METHOD
        ADD CONSTRAINT PK_METHOD PRIMARY KEY (METHOD_AC)
        """.format(schema)
    )
    cur.execute('CREATE INDEX I_METHOD$DBCODE ON {}.METHOD (DBCODE) NOLOGGING'.format(schema))
    oracledb.gather_stats(cur, schema, 'METHOD', cascade=True)
    oracledb.grant(cur, 'SELECT', schema, 'METHOD', 'INTERPRO_SELECT')
    cur.close()
    con.close()

    logging.info('loading member database signatures: complete')


def report_swiss_descriptions(dsn, schema, **kwargs):
    filepath = kwargs.get('output')

    if not isinstance(filepath, str) or not filepath:
        logging.critical("invalid 'output' argument")
        exit(1)

    dirname = os.path.dirname(os.path.abspath(filepath))
    if not os.path.dirname(dirname):
        os.makedirs(dirname)

    # Get descriptions for before update (METHOD2SWISS_DE was populated during protein update)
    con = oracledb.connect(dsn)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, M.DESCRIPTION
        FROM {0}.METHOD2SWISS_DE M 
        INNER JOIN {0}.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        """.format(schema)
    )

    old_entries = {}
    for acc, descr in cur:
        if acc in old_entries:
            old_entries[acc].append(descr)
        else:
            old_entries[acc] = [descr]

    # Get descriptions for after update
    cur.execute(
        """
        SELECT DISTINCT EM.ENTRY_AC, D.TEXT
          FROM {0}.METHOD2PROTEIN MP
        INNER JOIN {0}.METHOD M ON MP.METHOD_AC = M.METHOD_AC
        INNER JOIN {0}.DESC_VALUE D ON MP.DESC_ID = D.DESC_ID
        INNER JOIN {0}.ENTRY2METHOD EM ON M.METHOD_AC = EM.METHOD_AC
        WHERE MP.DBCODE = 'S'  
        """.format(schema)
    )

    new_entries = {}
    for acc, descr in cur:
        if acc in new_entries:
            new_entries[acc].append(descr)
        else:
            new_entries[acc] = [descr]

    cur.execute('SELECT ENTRY_AC, ENTRY_TYPE, CHECKED FROM INTERPRO_ANALYSIS.ENTRY')
    entry_info = {acc: (entry_type, is_checked) for acc, entry_type, is_checked in cur}

    cur.close()
    con.close()

    entries = {}
    for acc in old_entries:
        try:
            entry_type, is_checked = entry_info[acc]
        except KeyError:
            continue

        if acc in new_entries:
            new_entry_descrs = new_entries.pop(acc)
            entries[acc] = {
                'acc': acc,
                'type': entry_type,
                'checked': is_checked,
                'lost': [descr for descr in old_entries[acc] if descr not in new_entry_descrs],
                'gained': [descr for descr in new_entry_descrs if descr not in old_entries[acc]]
            }
        else:
            # All descriptions lost (as we lost the entry)
            entries[acc] = {
                'acc': acc,
                'type': entry_type,
                'checked': is_checked,
                'lost': old_entries[acc],
                'gained': []
            }

    for acc in new_entries:
        entry_type, is_checked = entry_info[acc]
        entries[acc] = {
            'acc': acc,
            'type': entry_type,
            'checked': is_checked,
            'lost': [],
            'gained': new_entries[acc]
        }

    with open(filepath, 'wt') as fh:
        fh.write('Entry\tType\tChecked\t#Gained\t#Lost\tGained\tLost\n')

        for entry in sorted(entries.values(), key=lambda e: (0 if e['type'] == 'F' else 1, e['type'], e['acc'])):
            fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                entry['acc'],
                entry['type'],
                entry['checked'],
                len(entry['gained']),
                len(entry['lost']),
                ' | '.join(entry['gained']),
                ' | '.join(entry['lost'])
            ))
