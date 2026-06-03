#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail query builder
#


class QueryBuilder:
    """
    Currently only configured to work with sqlite, will have to subclass once we add more database types
    """

    def __init__(self):
        self.selects = []
        self.from_table = None
        self.joins = []
        self.wheres = []
        self.group_bys = []
        self.order_by = None
        self.delete_from = None
        self.drop_if_exists = None
        self.subqueries = []
        self.params = []
        self.aliases = {}
        self.insert_into = None
        self.returning = None
        self.limit = None
        self.descending = False

    def INSERT_INTO(self, table_name, *columns):
        if columns:
            column_statement = f""" ({", ".join(map(str, columns))}) VALUES ({",".join(["?"] * len(columns))})"""
        else:
            column_statement = ""
        self.insert_into = f"INSERT INTO {table_name}{column_statement}"

        return self

    def RETURNING(self, column_name: str):
        self.returning = column_name
        return self

    def SELECT(self, *fields):
        self.selects.extend(fields)
        return self

    def SELECT_STATUS(self):
        status_case = """CASE
            WHEN EXISTS (SELECT 1 FROM Accepted s WHERE s.pose_id = R.pose_id) THEN 'accepted'
            WHEN EXISTS (SELECT 1 FROM Rejected s WHERE s.pose_id = R.pose_id) THEN 'rejected'
            WHEN EXISTS (SELECT 1 FROM Maybe s WHERE s.pose_id = R.pose_id) THEN 'maybe'
            ELSE ''
        END AS status"""
        self.selects.append(status_case)
        return self

    def FROM_BOOKMARK(self, bookmark, alias=None, db_alias=""):
        return self.FROM(f"({self.bookmark_query(bookmark, db_alias)})", alias)

    @staticmethod
    def bookmark_query(bookmark, db_alias=""):
        # bookmark names are validated to [a-z0-9_] by valid_bookmark_name() before reaching here
        if not all(c.isalnum() or c == "_" for c in bookmark):
            raise ValueError(f"Unsafe bookmark name: {bookmark!r}")
        if db_alias:
            db_alias += "."
        return (
            f"SELECT pose_id FROM {db_alias}filtered_poses "
            f"WHERE filter_id = "
            f"(SELECT filter_id FROM {db_alias}Filters WHERE name = '{bookmark}')"
        )

    def IN_BOOKMARK(self, bookmark):
        return self.WHERE(f"{self.aliased('Results')}.pose_id IN ({self.bookmark_query(bookmark)})")

    def FROM(self, table, alias=None, db_name=None):
        if db_name and alias:
            alias = db_name + "_" + alias
        if alias:
            alias = self._add_alias(table, alias)
        else:
            alias = self._add_alias(table, table)
        if db_name:
            table = db_name + "." + table
        self.from_table = (table, alias)
        return self

    def JOIN(self, table, alias=None, on: str = None, to=None, kind="INNER"):
        if on is None:
            raise ValueError(f"JOIN on table '{table}' requires an 'on' column")
        if alias:
            alias = self._add_alias(table, alias)
        if to:
            join_on = self.aliased(to) + "." + on + "=" + alias + "." + on
        else:
            join_on = (
                self.aliased(self.from_table[-1]) + "." + on + "=" + alias + "." + on
            )
        self.joins.append((kind, table, alias, join_on))
        return self

    def WHERE(self, condition, *params):
        self.wheres.append(condition)
        self.params.extend(params)
        return self

    def GROUP_BY(self, *fields):
        self.group_bys.extend(fields)
        return self

    def ORDER_BY(self, column):
        self.order_by = column
        return self

    def LIMIT(self, limit: int):
        self.limit = limit
        return self

    def DESC(self, descending):
        self.descending = descending
        return self

    def WITH_SUBQUERY(self, name, query, params=None):
        self.subqueries.append((name, query))
        if params:
            self.params.extend(params)
        return self

    def DELETE_FROM(self, table: str):
        self.delete_from = f"""DELETE FROM {table}"""
        return self

    def DROP_IF_EXISTS(self, table: str):
        self.drop_if_exists = f"DROP TABLE IF EXISTS {table}"
        return self

    def aliased(self, table: str):
        return self.aliases.get(table.lower(), table.lower())

    def _add_alias(self, table: str, alias: str) -> str:
        if alias.lower() in self.aliases.values():
            raise ValueError(f"Alias {alias} already in use for table {table}")
        if table.lower() in self.aliases:
            raise ValueError(
                f"Table {table} already has an alias: {self.aliases[table.lower()]}. Please use this instead."
            )
        self.aliases[table.lower()] = alias.lower()
        return alias

    def build(self, count=False):
        parts = []
        if self.insert_into:
            parts.append(f"""{self.insert_into}""")

        if self.subqueries:
            ctes = ", ".join(f"{name} AS ({query})" for name, query in self.subqueries)
            parts.append(f"WITH {ctes}")

        if self.selects:
            parts.extend(["SELECT", ", ".join(self.selects)])

        if self.delete_from:
            parts.append(self.delete_from)

        if self.drop_if_exists:
            parts.append(self.drop_if_exists)

        if self.from_table:
            table, alias = self.from_table
            parts.append("FROM {}{}".format(table, f" AS {alias}" if alias else ""))

        for kind, table, alias, join_on in self.joins:
            join_clause = f"{kind} JOIN {table}"
            if self.aliases.get(table.lower()):
                join_clause += f" AS {self.aliases[table.lower()]}"
            if join_on:
                join_clause += f" ON {join_on}"
            parts.append(join_clause)

        if self.wheres:
            parts.append("WHERE " + " AND ".join(self.wheres))

        if self.group_bys:
            parts.append("GROUP BY " + ", ".join(self.group_bys))

        if self.order_by:
            parts.append("ORDER BY " + self.order_by)

        if self.descending:
            parts.append("DESC")

        if self.limit is not None:
            parts.append(f"LIMIT {self.limit}")

        if self.returning:
            parts.append("RETURNING " + self.returning)

        sql = " ".join(parts)
        if count:
            if not self.selects:
                raise ValueError("build(count=True) is only valid for SELECT queries")
            return f"SELECT COUNT(*) FROM ({sql})", self.params
        return sql, self.params


class QueryBuilderSQLite(QueryBuilder):

    pass


class QueryBuilderDuck(QueryBuilder):
    """
    At some point I need to keep a wrapped SELECTS list so I can keep track,
    some times I might need to modify how it is wrapped
    """

    def GROUP_BY(self, *fields):
        self.group_bys.extend(fields)
        group_set = {f.lower() for f in self.group_bys}
        for index, item in enumerate(self.selects):
            if "ANY_VALUE" not in item and item.lower() not in group_set:
                self.selects[index] = f"ANY_VALUE({item})"
        return self
