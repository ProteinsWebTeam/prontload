#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cx_Oracle


def connect(dsn):
    return cx_Oracle.connect(dsn)


def drop_table(cur, ownname, tabname):
    try:
        cur.execute('DROP TABLE {}.{} CASCADE CONSTRAINTS'.format(ownname, tabname))
    except cx_Oracle.DatabaseError as e:
        """
        From cx_Oracle documentation:
        With cx_Oracle every exception object has exactly one argument in the `args` tuple.
        This argument is a `cx_Oracle._Error` object which has the following five read-only attributes.

        (http://cx-oracle.readthedocs.io/en/latest/module.html?highlight=DatabaseError)
        """
        _error = e.args[0]
        if _error.code == 942:
            # ORA-00942 (table or view does not exist)
            # That's fine since we are going to create the table
            pass
        else:
            # Something else: raise the issue
            raise e


def gather_stats(cur, ownname, tabname, cascade=True):
    if cascade:
        cur.execute(
            """
            begin
                dbms_stats.gather_table_stats(:1, :2, cascade=>TRUE);
            end;
            """,
            (ownname, tabname)
        )
    else:
        cur.callproc('DBMS_STATS.GATHER_TABLE_STATS', (ownname, tabname))


def create_synonym(cur, src, dst, obj):
    cur.execute('CREATE OR REPLACE SYNONYM {0}.{2} FOR {1}.{2}'.format(dst, src, obj))


def list_tables(cur, ownname):
    cur.execute(
        """
        SELECT
          table_name
        FROM
          dba_tables
        WHERE
          owner=:1
        ORDER BY
          table_name
        """,
        (ownname,)
    )

    return [row[0] for row in cur]


def grant(cur, privilege, schema, tabname, grantee):
    cur.execute('GRANT {} ON {}.{} TO {}'.format(privilege, schema, tabname, grantee))
