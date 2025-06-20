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
        self.subquery = {}
        self.params = []
        self.aliases = {}
        self.insert_into = None
        self.limit = None

    def INSERT_INTO(self, table_name):
        self.insert_into = f"INSERT INTO {table_name}"
        return self

    def SELECT(self, *fields):
        self.selects.extend(fields)
        return self

    def SELECT_STATUS(self):
        status_case = """CASE
            WHEN EXISTS (SELECT 1 FROM Accepted s WHERE s.pose_id = R.pose_ID) THEN 'accepted'
            WHEN EXISTS (SELECT 1 FROM Rejected s WHERE s.pose_id = R.pose_ID) THEN 'rejected'
            WHEN EXISTS (SELECT 1 FROM Maybe s WHERE s.pose_id = R.pose_ID) THEN 'maybe'
            ELSE 'not evaluated'
        END AS status"""
        self.selects.append(status_case)
        return self

    def FROM_BOOKMARK(self, bookmark, alias=None, db_alias=""):
        bookmark_query = f"({self.bookmark_query(bookmark, db_alias)})"
        return self.FROM(bookmark_query, alias)

    @staticmethod
    def bookmark_query(bookmark, db_alias=""):
        if db_alias:
            db_alias += "."
        query = f"""SELECT Pose_id FROM {db_alias}filtered_poses 
        WHERE filter_id = 
            (SELECT filter_id FROM {db_alias}Filters 
            WHERE name = '{bookmark}')"""

        return query

    def IN_BOOKMARK(self, bookmark):
        """
        should make string that adds to
        <results-alias> in (<bookmark_query>)
        """
        query = f"""{self.aliased("Results")}.pose_id IN ({self.bookmark_query(bookmark)})"""
        return self.WHERE(query)

    def FROM(self, table, alias=None):
        if alias:
            alias = self._add_alias(table, alias)
        self.from_table = (table, alias)
        return self

    def JOIN(self, table, alias=None, on=None, to=None, kind="INNER"):
        if alias:
            alias = self._add_alias(table, alias)
        if to:
            join_on = self.aliased(to) + "." + on + "=" + alias + "." + on
        else:
            # use last added select table
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
        self.limit = str(limit)

    def WITH_SUBQUERY(self, name, query, params=None):
        # This could in theory accept a QueryBuilder object instead
        self.subquery["name"] = name
        self.subquery["query"] = query
        if params:
            self.params.extend(params)
        return self

    def aliased(self, table: str):
        if table.lower() in self.aliases.keys():
            return self.aliases[table.lower()]
        else:
            return table.lower()

    def _add_alias(self, table: str, alias: str) -> str:
        if alias.lower() in self.aliases.values():
            raise Exception(f"Alias {alias} already in use for table {table}")
        if table.lower() in self.aliases.keys():
            raise Exception(
                f"WARNING: Table {table} already has an alias: {self.aliases[table]}. Please use this instead."
            )
        self.aliases[table.lower()] = alias.lower()
        return alias

    def build(self, count=False):
        parts = []
        if self.insert_into:
            parts.append(self.insert_into)

        if self.subquery:
            parts.append(
                f"""WITH {self.subquery["name"]} AS ({self.subquery["query"]})"""
            )
        parts.extend(["SELECT", ", ".join(self.selects)])

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

        if self.limit:
            parts.append("LIMIT " + self.limit)

        if count:
            return f"SELECT COUNT(*) FROM ({' '.join(parts)})", self.params
        else:
            return " ".join(parts), self.params
