#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail database schema and supporting functionality
#


from dataclasses import dataclass, field
from typing import Union

NUMERIC_TYPES = {"INTEGER", "FLOAT"}


@dataclass
class Column:
    sql_type: str
    description: str
    primary_key: bool = False
    autoincrement: bool = False
    nullable: bool = True
    unique: bool = False
    on_conflict_ignore: bool = False  # SQLite only: ON CONFLICT IGNORE on UNIQUE
    foreign_key: Union[str, None] = None  # "ReferencedTable.column"
    default: Union[str, None] = None  # raw SQL default expression


@dataclass
class TableSchema:
    columns: dict[str, Column]
    unique_together: list[list[str]] = field(default_factory=list)
    # Each inner list is one multi-column UNIQUE constraint.
    # SQLite emits ON CONFLICT IGNORE on these automatically.
    sqlite_indices: list[list[str]] = field(default_factory=list)
    # SQLite index column groups. DuckDB ignores these (auto-optimises).
    temporary: bool = False
    name: str = ""


# ---------------------------------------------------------------------------
# Dialect type maps
# ---------------------------------------------------------------------------

SQLITE_TYPES: dict[str, str] = {
    "INTEGER": "INTEGER",
    "FLOAT": "FLOAT",
    "VARCHAR": "VARCHAR",
    "BLOB": "BLOB",
    "DATETIME": "DATETIME",
    "JSON_OR_VARCHAR": "VARCHAR",
}

DUCKDB_TYPES: dict[str, str] = {
    "INTEGER": "INTEGER",
    "FLOAT": "FLOAT",
    "VARCHAR": "VARCHAR",
    "BLOB": "BLOB",
    "DATETIME": "TIMESTAMP",
    "JSON_OR_VARCHAR": "JSON",
}

# ---------------------------------------------------------------------------
# Table schemas
# ---------------------------------------------------------------------------

LIGANDS_SCHEMA = TableSchema(
    name="Ligands",
    columns={
        "ligand_id": Column(
            "INTEGER",
            "unique ligand identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "ligname": Column(
            "VARCHAR",
            "ligand name",
            nullable=False,
            unique=True,
            on_conflict_ignore=True,
        ),
        "ligand_smile": Column("VARCHAR", "SMILES string"),
        "rdmol": Column("BLOB", "RDKit molecule binary"),
    },
)

RESULTS_SCHEMA = TableSchema(
    name="Results",
    columns={
        "pose_id": Column(
            "INTEGER",
            "unique pose identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "ligand_id": Column(
            "INTEGER", "owning ligand", foreign_key="Ligands.ligand_id"
        ),
        "receptor": Column("VARCHAR", "receptor name"),
        "pose_rank": Column("INTEGER", "rank of ligand pose"),
        "run_number": Column("INTEGER", "run number for ligand pose"),
        "docking_score": Column("FLOAT", "docking score / energy"),
        "leff": Column("FLOAT", "ligand efficiency"),
        "delta": Column("FLOAT", "delta energy from best pose"),
        "cluster_rmsd": Column("FLOAT", "RMSD within cluster"),
        "cluster_size": Column("INTEGER", "cluster size"),
        "reference_rmsd": Column("FLOAT", "RMSD to reference pose"),
        "energies_inter": Column("FLOAT", "intermolecular energy"),
        "energies_vdw": Column("FLOAT", "van der Waals energy"),
        "energies_electro": Column("FLOAT", "electrostatic energy"),
        "energies_flexLig": Column("FLOAT", "flexible ligand energy"),
        "energies_flexLR": Column("FLOAT", "flexible receptor energy"),
        "energies_intra": Column("FLOAT", "intramolecular energy"),
        "energies_torsional": Column("FLOAT", "torsional energy"),
        "unbound_energy": Column("FLOAT", "unbound state energy"),
        "num_interactions": Column("INTEGER", "number of interactions"),
        "num_hb": Column("INTEGER", "hydrogen bonds"),
        "pose_coordinates": Column("VARCHAR", "pose atom coordinates"),
        "flexible_res_coordinates": Column("VARCHAR", "flexible residue coordinates"),
    },
    sqlite_indices=[
        ["docking_score", "leff"],
        ["ligand_id"],
    ],
)

RECEPTORS_SCHEMA = TableSchema(
    name="Receptors",
    columns={
        "receptor_id": Column(
            "INTEGER",
            "unique receptor identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "recname": Column("VARCHAR", "receptor name"),
        "box_dim": Column("VARCHAR", "docking box dimensions"),
        "box_center": Column("VARCHAR", "docking box center"),
        "grid_spacing": Column("FLOAT", "grid spacing"),
        "flexible_residues": Column("VARCHAR", "flexible residue names"),
        "flexres_atomnames": Column("VARCHAR", "flexible residue atom names"),
        "receptor_object": Column("BLOB", "receptor binary object"),
        "polymer": Column("JSON_OR_VARCHAR", "polymer data"),
    },
)

DB_PROPERTIES_SCHEMA = TableSchema(
    name="db_properties",
    columns={
        "db_write_session": Column(
            "INTEGER",
            "write session identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "docking_mode": Column("VARCHAR", "docking mode used"),
        "number_of_poses": Column("INTEGER", "number of poses in session"),
    },
)

INTERACTION_INDICES_SCHEMA = TableSchema(
    name="Interaction_indices",
    columns={
        "interaction_id": Column(
            "INTEGER",
            "unique interaction type identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "interaction_type": Column("VARCHAR", "type of interaction"),
        "rec_chain": Column("VARCHAR", "receptor chain"),
        "rec_resname": Column("VARCHAR", "receptor residue name"),
        "rec_resid": Column("VARCHAR", "receptor residue id"),
        "rec_atom": Column("VARCHAR", "receptor atom name"),
        "rec_atomid": Column("VARCHAR", "receptor atom id"),
    },
    unique_together=[
        [
            "interaction_type",
            "rec_chain",
            "rec_resname",
            "rec_resid",
            "rec_atom",
            "rec_atomid",
        ]
    ],
    sqlite_indices=[
        [
            "interaction_type",
            "rec_chain",
            "rec_resname",
            "rec_resid",
            "rec_atom",
            "rec_atomid",
        ]
    ],
)

INTERACTIONS_SCHEMA = TableSchema(
    name="Interactions",
    columns={
        "interaction_pose_id": Column(
            "INTEGER",
            "unique interaction-pose record",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "pose_id": Column("INTEGER", "pose", foreign_key="Results.pose_id"),
        "interaction_id": Column(
            "INTEGER",
            "interaction index",
            foreign_key="Interaction_indices.interaction_id",
        ),
    },
    sqlite_indices=[["pose_id", "interaction_id"]],
)

FILTERS_SCHEMA = TableSchema(
    name="Filters",
    columns={
        "filter_id": Column(
            "INTEGER",
            "unique filter/bookmark identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "name": Column("VARCHAR", "bookmark name"),
        "query": Column("VARCHAR", "SQL query string"),
        "filters": Column("VARCHAR", "filter parameters as string"),
        "filter_window": Column("VARCHAR", "filter window parameters"),
    },
)

FILTERED_POSES_SCHEMA = TableSchema(
    name="Filtered_poses",
    columns={
        "filter_id": Column(
            "INTEGER", "bookmark reference", foreign_key="Filters.filter_id"
        ),
        "pose_id": Column("INTEGER", "pose reference", foreign_key="Results.pose_id"),
    },
)

CLUSTERS_SCHEMA = TableSchema(
    name="Clusters",
    columns={
        "cluster_id": Column(
            "INTEGER",
            "unique clustering session identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "name": Column("VARCHAR", "cluster name"),
        "description": Column("VARCHAR", "cluster description"),
        "cluster_window": Column("VARCHAR", "clustering parameters"),
        "num_clusters": Column("INTEGER", "number of clusters"),
    },
)

CLUSTER_GROUPS_SCHEMA = TableSchema(
    name="Cluster_groups",
    columns={
        "cluster_id": Column(
            "INTEGER", "clustering session", foreign_key="Clusters.cluster_id"
        ),
        "cluster_group": Column("INTEGER", "cluster group number"),
        "representative": Column(
            "INTEGER", "representative pose", foreign_key="Results.pose_id"
        ),
    },
)

POSE_CLUSTERS_SCHEMA = TableSchema(
    name="Pose_clusters",
    columns={
        "cluster_id": Column(
            "INTEGER", "clustering session", foreign_key="Clusters.cluster_id"
        ),
        "cluster_group": Column("INTEGER", "cluster group number"),
        "pose_id": Column("INTEGER", "pose", foreign_key="Results.pose_id"),
    },
)

MERGED_TABLES_SCHEMA = TableSchema(
    name="Merged_tables",
    columns={
        "merge_id": Column(
            "INTEGER",
            "unique merge session identifier",
            primary_key=True,
            autoincrement=True,
            nullable=False,
        ),
        "dbfile": Column("VARCHAR", "path to merged database file"),
        "merge_start": Column(
            "DATETIME", "merge start timestamp", default="CURRENT_TIMESTAMP"
        ),
    },
)

PK_CONVERSIONS_SCHEMA = TableSchema(
    name="PK_conversions",
    columns={
        "merge_id": Column(
            "INTEGER", "merge session", foreign_key="merged_tables.merge_id"
        ),
        "table_name": Column("VARCHAR", "table with remapped primary keys"),
        "original_PK": Column("INTEGER", "original primary key value"),
        "merged_PK": Column("INTEGER", "new primary key after merge"),
    },
    sqlite_indices=[["merge_id", "original_PK"]],
)

# Accepted, Maybe, Rejected share this structure; pass the name to build_create_table.
STATUS_TABLE_SCHEMA = TableSchema(
    columns={
        "pose_id": Column(
            "INTEGER", "pose", unique=True, foreign_key="Results.pose_id"
        ),
    }
)

CANDIDATES_SUBQ = "(SELECT pose_id FROM Accepted UNION SELECT pose_id FROM Maybe)"
CANDIDATES_NAME = "candidates"
LIGANDS_ONLY_COLS = set(LIGANDS_SCHEMA.columns) - set(RESULTS_SCHEMA.columns)

# Table of comments from e.g., visual inspection
POSE_COMMENTS_SCHEMA = TableSchema(
    name="Pose_comments",
    columns={
        "pose_id": Column(
            "INTEGER", "pose", primary_key=True, foreign_key="Results.pose_id"
        ),
        "comment": Column("VARCHAR", "user comment on pose"),
    }
)

TABLE_SCHEMAS: dict[str, TableSchema] = {
    s.name.lower(): s
    for s in [
        LIGANDS_SCHEMA,
        RESULTS_SCHEMA,
        RECEPTORS_SCHEMA,
        DB_PROPERTIES_SCHEMA,
        INTERACTION_INDICES_SCHEMA,
        INTERACTIONS_SCHEMA,
        FILTERS_SCHEMA,
        FILTERED_POSES_SCHEMA,
        CLUSTERS_SCHEMA,
        CLUSTER_GROUPS_SCHEMA,
        POSE_CLUSTERS_SCHEMA,
        MERGED_TABLES_SCHEMA,
        PK_CONVERSIONS_SCHEMA,
        POSE_COMMENTS_SCHEMA,
    ]
}
# Status tables share a schema — register them separately
for _status in ("accepted", "maybe", "rejected"):
    TABLE_SCHEMAS[_status] = STATUS_TABLE_SCHEMA

ALL_TABLE_NAMES = frozenset(TABLE_SCHEMAS)

# ---------------------------------------------------------------------------
# Derived schemas for outfields and order_results
# ---------------------------------------------------------------------------

# All columns a user can request in output — internal keys excluded.
OUTFIELD_SCHEMA: dict[str, Column] = {
    name: col
    for schema in [
        RESULTS_SCHEMA,
        LIGANDS_SCHEMA,
        INTERACTION_INDICES_SCHEMA,
    ]
    for name, col in schema.columns.items()
    if not col.primary_key and not col.foreign_key
}

# Same columns grouped by source table — used for export and column-to-table resolution.
OUTFIELD_BY_TABLE: dict[str, dict[str, Column]] = {
    "Results": {
        n: c
        for n, c in RESULTS_SCHEMA.columns.items()
        if not c.primary_key and not c.foreign_key
    },
    "Ligands": {
        n: c
        for n, c in LIGANDS_SCHEMA.columns.items()
        if not c.primary_key and not c.foreign_key
    },
    "Interaction_indices": {
        n: c
        for n, c in INTERACTION_INDICES_SCHEMA.columns.items()
        if not c.primary_key and not c.foreign_key
    },
}

# Only numeric columns make sense as an ORDER BY key.
ORDER_RESULT_SCHEMA: dict[str, Column] = {
    name: col
    for name, col in OUTFIELD_SCHEMA.items()
    if col.sql_type in ("FLOAT", "INTEGER", "VARCHAR")
}

# ---------------------------------------------------------------------------
# DDL builder
# ---------------------------------------------------------------------------


def build_create_table(table_name: str, schema: TableSchema, dialect: str) -> list[str]:
    """Return SQL statements that create the table for the given dialect.

    May return more than one statement (e.g. DuckDB emits a CREATE SEQUENCE
    before the CREATE TABLE for each autoincrement column).

    Args:
        table_name: SQL table name.
        schema: TableSchema instance describing the table.
        dialect: "sqlite" or "duckdb".

    Returns:
        list[str]: Ordered SQL statements to execute.
    """
    type_map = SQLITE_TYPES if dialect == "sqlite" else DUCKDB_TYPES
    statements: list[str] = []
    col_defs: list[str] = []
    fk_constraints: list[str] = []

    for col_name, col in schema.columns.items():
        type_str = type_map[col.sql_type]
        parts = [col_name, type_str]

        if col.autoincrement:
            if dialect == "sqlite":
                parts.append("PRIMARY KEY AUTOINCREMENT")
            else:
                seq = f"seq_{table_name.lower()}_{col_name.lower()}"
                statements.append(f"CREATE SEQUENCE IF NOT EXISTS {seq} START 1")
                parts.append(f"DEFAULT nextval('{seq}') PRIMARY KEY")
        elif col.primary_key:
            parts.append("PRIMARY KEY")

        if not col.nullable and not col.primary_key and not col.autoincrement:
            parts.append("NOT NULL")

        if col.unique:
            if col.on_conflict_ignore and dialect == "sqlite":
                parts.append("UNIQUE ON CONFLICT IGNORE")
            else:
                parts.append("UNIQUE")

        if col.default:
            parts.append(f"DEFAULT {col.default}")

        if col.foreign_key:
            ref_table, ref_col = col.foreign_key.split(".")
            if dialect == "sqlite":
                fk_constraints.append(
                    f"FOREIGN KEY ({col_name}) REFERENCES {ref_table}({ref_col})"
                )
            else:
                parts.append(f"REFERENCES {ref_table}({ref_col})")

        col_defs.append(" ".join(parts))

    for unique_cols in schema.unique_together:
        constraint = f"UNIQUE ({', '.join(unique_cols)})"
        if dialect == "sqlite":
            constraint += " ON CONFLICT IGNORE"
        col_defs.append(constraint)

    all_defs = col_defs + fk_constraints
    prefix = "CREATE TEMP TABLE" if schema.temporary else "CREATE TABLE"
    statements.append(f"{prefix} IF NOT EXISTS {table_name} ({', '.join(all_defs)})")
    return statements


def build_create_indices(table_name: str, schema: TableSchema) -> list[str]:
    """Return CREATE INDEX statements for a table (SQLite only).

    Args:
        table_name: SQL table name.
        schema: TableSchema instance.

    Returns:
        list[str]: CREATE INDEX statements, empty list if schema has no indices.
    """
    statements = []
    for i, index_cols in enumerate(schema.sqlite_indices):
        index_name = (
            f"ak_{table_name.lower()}_{i}" if i > 0 else f"ak_{table_name.lower()}"
        )
        cols = ", ".join(index_cols)
        statements.append(
            f"CREATE INDEX IF NOT EXISTS {index_name} ON {table_name}({cols})"
        )
    return statements
