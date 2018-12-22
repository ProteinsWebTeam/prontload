#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cx_Oracle


BULK_INSERT_SIZE = 100000


class Connection(object):
    def __init__(self, dsn):
        self._dsn = dsn
        self._db = None
        self.open()

    def open(self):
        self.close()
        self._db = cx_Oracle.connect(self._dsn, encoding='utf-8',
                                     nencoding='utf-8')

    def close(self):
        if self._db is not None:
            self._db.close()
            self._db = None

    def commit(self):
        if self._db is not None:
            self._db.commit()

    def get(self, query, *args, **kwargs):
        cur = self._cursor()
        try:
            self._execute(cur, query, args, kwargs)
            #column_names = [col[0] for col in cur.description]
            for row in cur:
                # yield dict(zip(column_names, row))
                yield row
        finally:
            cur.close()

    def exec(self, procedure, *args):
        cur = self._cursor()
        try:
            cur.callproc(procedure, args)
        finally:
            cur.close()

    def execute(self, query, *args, **kwargs):
        cur = self._cursor()
        try:
            self._execute(cur, query, args, kwargs)
            return cur.rowcount
        finally:
            cur.close()

    def executemany(self, query, args):
        cur = self._cursor()
        try:
            self._executemany(cur, query, args)
            return cur.rowcount
        finally:
            cur.close()

    def get_tables(self, ownname):
        query = """
            SELECT
              table_name
            FROM
              dba_tables
            WHERE
              UPPER(owner) = :1
            ORDER BY
              table_name
        """
        return map(lambda r: r[0], self.get(query, ownname))

    def drop_table(self, ownname, tabname, forgive_busy=False):
        try:
            self.execute("DROP TABLE {}.{} "
                         "CASCADE CONSTRAINTS".format(ownname, tabname))
        except cx_Oracle.DatabaseError as e:
            """
            From cx_Oracle documentation:
            With cx_Oracle every exception object has exactly one argument
                in the `args` tuple.
            This argument is a `cx_Oracle._Error` object
                which has the following five read-only attributes.

            (http://cx-oracle.readthedocs.io/en/latest/module.html?highlight=DatabaseError)
            """
            _error = e.args[0]
            if _error.code == 942:
                # ORA-00942 (table or view does not exist)
                # That's fine since we are going to create the table
                return False
            elif _error.code == 54 and forgive_busy:
                """
                ORA-00054
                resource busy and acquire with NOWAIT specified
                    or timeout expired
                """
                return False
            else:
                # Something else: raise the issue
                raise e
        else:
            return True

    def optimise_table(self, ownname, tabname, cascade=False):
        if cascade:
            self.execute(
                """
                begin
                    dbms_stats.gather_table_stats(:1, :2, cascade=>TRUE);
                end;
                """,
                ownname, tabname
            )
        else:
            self.exec("DBMS_STATS.GATHER_TABLE_STATS", (ownname, tabname))

    def grant(self, privilege, schema, tabname, grantee):
        self.execute("GRANT {} "
                     "ON {}.{} "
                     "TO {}".format(privilege, schema, tabname, grantee))

    def _cursor(self):
        return self._db.cursor()

    @staticmethod
    def _execute(cursor, query, args, kwargs):
        return cursor.execute(query, args or kwargs)

    @staticmethod
    def _executemany(cursor, query, args):
        return cursor.executemany(query, args)

    def __del__(self):
        self.close()
