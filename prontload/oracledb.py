from concurrent.futures import as_completed, ThreadPoolExecutor

import cx_Oracle

from . import get_logger

logger = get_logger()


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
            SELECT TABLE_NAME 
            FROM ALL_TABLES
            WHERE UPPER(OWNER) = UPPER(:1)
            ORDER BY TABLE_NAME
        """
        return [row[0] for row in self.get(query, ownname)]

    def get_indexes(self, ownname):
        query = """
            SELECT 
              I.INDEX_NAME, I.TABLESPACE_NAME, I.UNIQUENESS, I.TABLE_NAME, 
              I.LOGGING, IC.COLUMN_NAME, IC.DESCEND
            FROM ALL_INDEXES I
            INNER JOIN ALL_IND_COLUMNS IC 
              ON I.OWNER = IC.INDEX_OWNER 
              AND I.INDEX_NAME = IC.INDEX_NAME 
              AND I.TABLE_NAME = IC.TABLE_NAME
            WHERE UPPER(I.OWNER) = UPPER(:1)
            ORDER BY I.INDEX_NAME, IC.COLUMN_POSITION
        
        """
        indexes = {}
        for row in self.get(query, ownname):
            name = row[0]
            if name in indexes:
                idx = indexes[name]
            else:
                idx = indexes[name] = {
                    "name": name,
                    "tablespace": row[1],
                    "is_unique": row[2] == "UNIQUE",
                    "table": row[3],
                    "logging": row[4],
                    "columns": []
                }

            idx["columns"].append({
                "name": row[5],
                "order":  row[6]
            })

        return list(indexes.values())

    def get_constraints(self, ownname):
        query = """
            SELECT 
              C.CONSTRAINT_NAME, C.CONSTRAINT_TYPE, C.TABLE_NAME,
              C.R_CONSTRAINT_NAME, CC.COLUMN_NAME
            FROM ALL_CONSTRAINTS C
            INNER JOIN ALL_CONS_COLUMNS CC 
              ON C.OWNER = CC.OWNER 
              AND C.CONSTRAINT_NAME = CC.CONSTRAINT_NAME 
              AND C.TABLE_NAME = CC.TABLE_NAME
            WHERE UPPER(C.OWNER) = UPPER(:1)
            ORDER BY C.CONSTRAINT_NAME, CC.POSITION
        """

        constraints = {}
        for row in self.get(query, ownname):
            name = row[0]
            _type = row[1]

            if _type == 'P':
                _type = "PRIMARY KEY"
            elif _type == 'U':
                _type = "UNIQUE"
            else:
                # Only support primary/unique keys
                continue

            if name in constraints:
                const = constraints[name]
            else:
                const = constraints[name] = {
                    "name": name,
                    "table": row[2],
                    "type": _type,
                    "reference": row[3],
                    "columns": []
                }

            const["columns"].append(row[4])

        return list(constraints.values())

    def get_grants(self, owname):
        query = """
            SELECT TABLE_NAME, PRIVILEGE, GRANTEE 
            FROM ALL_TAB_PRIVS 
            WHERE TABLE_SCHEMA = :1
        """

        tables = {}
        for table, privilege, grantee in self.get(query, owname):
            if table in tables:
                t = tables[table]
            else:
                t = tables[table] = []

            t.append({
                "privilege": privilege,
                "grantee": grantee
            })

        return tables

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


def clear_schema(dsn, schema):
    con = Connection(dsn)
    for table in con.get_tables(schema):
        con.drop_table(schema, table)


def copy_tables(dsn, schema_src, schema_dst):
    clear_schema(dsn, schema_dst)

    con = Connection(dsn)
    jobs = {}
    for table in con.get_tables(schema_src):
        jobs[table] = {
            "indexes": [],
            "constraints": [],
            "grants": []
        }

    for table, grants in con.get_grants(schema_src).items():
        if table in jobs:
            jobs[table]["grants"] = grants

    for idx in con.get_indexes(schema_src):
        table = idx["table"]
        if table in jobs:
            jobs[table]["indexes"].append(idx)

    for const in con.get_constraints(schema_src):
        table = const["table"]
        if table in jobs:
            jobs[table]["constraints"].append(const)

    num_errors = 0
    with ThreadPoolExecutor() as executor:
        fs = {}

        for table, t in jobs.items():
            f = executor.submit(rebuild_table, dsn, table, schema_src,
                                schema_dst, t["constraints"], t["indexes"],
                                t["grants"])
            fs[f] = table

        for f in as_completed(fs):
            table = fs[f]

            try:
                f.result()
            except Exception as exc:
                logger.error("{} failed: {}".format(table, exc))
                num_errors += 1
            else:
                logger.info("{} done".format(table))

    return num_errors > 0


def rebuild_table(dsn, name, src, dst, constraints, indexes, grants):
    con = Connection(dsn)
    logger.debug("{}: creating table".format(name))
    con.drop_table(dst, name)
    con.execute(
        """
        CREATE TABLE {2}.{0} 
        NOLOGGING 
        AS 
        SELECT * FROM {1}.{0}
        """.format(name, src, dst)
    )

    for grant in grants:
        con.grant(grant["privilege"], dst, name, grant["grantee"])

    for const in constraints:
        logger.debug("{}: adding constraint {}".format(name, const["name"]))
        con.execute(
            """
            ALTER TABLE {0}.{1}
            ADD CONSTRAINT {2} {3} ({4})
            """.format(dst, name, const["name"], const["type"], ", ".join(const["columns"]))
        )

    for idx in indexes:
        columns = ["{name} {order}".format(**col) for col in idx["columns"]]

        try:
            logger.debug("{}: creating index {}".format(name, idx["name"]))
            con.execute(
                """
                  CREATE INDEX {0}.{1}
                  ON {0}.{2}({3}) NOLOGGING
                """.format(dst, idx["name"], name, ", ".join(columns))
            )
        except cx_Oracle.DatabaseError as exc:
            _error = exc.args[0]
            if _error.code == 955:
                """
                ORA-00955: name is already used by an existing object
                -> index was created when adding constraint
                """
                logger.debug("{}: skipping index {}".format(name, idx["name"]))
            else:
                raise exc
