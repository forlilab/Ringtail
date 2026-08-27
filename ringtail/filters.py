#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail filters
#
# Filtering is Ringtail's core capability, so the types describing a filter live here
# rather than alongside the general option objects. A whole filter specification is a
# `Filters` tree of AND/OR nodes over `Filter` leaves; a plain flat filter set is just a
# tree with one leaf, so both kinds render through one path in the storage manager.
#

from pathlib import Path
from typing import Callable, Union
from dataclasses import dataclass, asdict, field

from .exceptions import OptionError
from .logutils import get_logger
from .schema import (
    HB_COUNT_COLUMN,
    INTERACTION_TYPE_LETTERS,
    NUMERIC_FILTER_SCHEMA,
)

logger = get_logger(__name__)


@dataclass
class Filter:
    """One leaf of a filter expression: a set of criteria that are combined with AND.

    Internal — you don't construct these directly. :class:`Filters` builds them, either
    from flat keyword arguments or from a nested expression dict. They surface when
    reading a filter specification back, via ``Filters.leaves()`` and ``Filters.leaf``.
    """

    eworst: float = None
    ebest: float = None
    lebest: float = None
    leworst: float = None
    score_percentile: float = None
    le_percentile: float = None
    vdw_interactions: list = field(default_factory=list)
    hb_interactions: list = field(default_factory=list)
    reactive_interactions: list = field(default_factory=list)
    hb_count: int = None
    react_any: bool = None
    max_miss: int = 0
    ligand_name: str = None
    ligand_name_file: str = None
    ligand_operator: str = None
    ligand_substruct: str = None
    ligand_substruct_pos: list = None
    ligand_max_atoms: int = None
    ligand_min_molweight: float = None
    ligand_max_molweight: float = None

    def __post_init__(self):
        if self.vdw_interactions is None:
            self.vdw_interactions = []
        if self.hb_interactions is None:
            self.hb_interactions = []
        if self.reactive_interactions is None:
            self.reactive_interactions = []
        if self.max_miss is None:
            # "not given" means "miss nothing"; a None-check rather than int(bool(...)),
            # which would collapse max_miss=2 to 1
            self.max_miss = 0
        self.checks()

    def asdict(self):
        return asdict(self)

    def checks(self):
        """Ensures all values are internally consistent and valid. Runs once after all values are
        set initially"""
        if self.eworst is not None and self.score_percentile is not None:
            raise OptionError(
                "Cannot use 'eworst' together with 'score_percentile': they are two "
                "different cutoffs on docking score. Set one or the other."
            )

        if self.leworst is not None and self.le_percentile is not None:
            raise OptionError(
                "Cannot use 'leworst' together with 'le_percentile': they are two "
                "different cutoffs on ligand efficiency. Set one or the other."
            )

        if self.score_percentile is not None and (
            self.score_percentile < 0 or self.score_percentile > 100
        ):
            raise OptionError(
                f"Given 'score_percentile' {self.score_percentile} not allowed. Should be within percentile range of 0-100."
            )

        if self.le_percentile is not None and (
            self.le_percentile < 0 or self.le_percentile > 100
        ):
            raise OptionError(
                f"Given 'le_percentile' {self.le_percentile} not allowed. Should be within percentile range of 0-100."
            )

        if self.ligand_operator not in ["OR", "AND"] and (
            self.ligand_substruct or self.ligand_substruct_pos
        ):
            logger.debug(f"'ligand_operator' set to default 'OR'.")
            self.ligand_operator = "OR"

        if self.max_miss < 0:
            raise OptionError("'max_miss' must be greater than or equal to 0.")

        # normalize ligand names to one flat list so the count limit is meaningful
        if self.ligand_name:
            from .util import iterate_nested

            self.ligand_name = [
                name
                for item in iterate_nested(self.ligand_name)
                for name in item.split(",")
            ]
            if len(self.ligand_name) > 50:
                raise OptionError(
                    "The number of provided ligand names is too large, please prepare as a csv and use 'ligand_name_file' instead."
                )

        if (
            self.ligand_name_file
            and Path(self.ligand_name_file).suffix.lower() != ".csv"
        ):
            raise OptionError(
                f"The file of ligand names needs to be a csv file, cannot proceed with {Path(self.ligand_name_file).suffix.lower()}."
            )

    def to_criteria(self, percentile_cutoff: Callable = None) -> dict:
        """Group these filters into the structured criteria a query builder renders.

        Pure option-domain work: maps each filter onto the column it constrains (per
        ``schema.NUMERIC_FILTER_SCHEMA``), buckets criteria by kind, and normalizes values.
        Emits no SQL — the storage manager turns the ``(column, operator, value)`` tuples
        into dialect-specific predicates.

        Percentile filters are the exception: resolving a percentile to a cutoff requires
        reading the data, so the caller injects a resolver rather than this class holding a
        database handle.

        Args:
            percentile_cutoff (Callable[[float, str], float], optional): resolves
                (percentile, column) to a concrete cutoff. Required only when
                'score_percentile' or 'le_percentile' is set.

        Returns:
            dict: any of ``numeric`` (list of ``(column, operator, value)``), ``interactions``
            (list of ``[type_letter, chain, resname, resid, atom, wanted]``), ``max_miss``
            (present only alongside ``interactions``), and ``ligand`` (in-memory ligand
            filters, applied by RDKit rather than SQL). Empty buckets are omitted, so an
            empty dict means no filters were set.
        """
        numeric = []
        interactions = []
        ligand = {}
        max_miss = 0
        interaction_keys = self.get_filter_keys("interaction")
        ligand_keys = self.get_filter_keys("ligand")

        for key, value in self.asdict().items():
            if value is None:
                continue

            if key in NUMERIC_FILTER_SCHEMA:
                column, operator = NUMERIC_FILTER_SCHEMA[key]
                if key.endswith("_percentile"):
                    if percentile_cutoff is None:
                        raise OptionError(
                            f"Filtering on '{key}' needs a percentile resolver; "
                            "pass 'percentile_cutoff' to to_criteria()."
                        )
                    value = percentile_cutoff(value, column)
                numeric.append((column, operator, value))

            elif key == "hb_count":
                # positive means "at least this many"; negative means "no more than",
                # and 0 lands here so that "no hydrogen bonds" stays expressible
                if value > 0:
                    numeric.append((HB_COUNT_COLUMN, ">=", value))
                else:
                    numeric.append((HB_COUNT_COLUMN, "<=", -value))

            elif key in interaction_keys:
                for interact in value:
                    # interact is ["chain:res:resno:resatom", wanted(bool)]
                    interaction_string = (
                        INTERACTION_TYPE_LETTERS[key] + ":" + interact[0]
                    )
                    interactions.append(interaction_string.split(":") + [interact[1]])

            elif key == "react_any" and value:
                interactions.append(["R", "", "", "", "", True])

            elif key == "max_miss":
                max_miss = value

            elif key in ligand_keys or key == "ligand_name_file":
                if key == "ligand_substruct_pos":
                    for f in value:
                        f[1] = int(f[1])
                        for index in range(2, 6):
                            f[index] = float(f[index])
                ligand[key] = value

        criteria = {}
        if numeric:
            criteria["numeric"] = numeric
        if interactions:
            criteria["interactions"] = interactions
            criteria["max_miss"] = max_miss
        if ligand:
            criteria["ligand"] = ligand
        return criteria

    # criteria RDKit has to evaluate in memory — they can't become SQL predicates
    RDKIT_CRITERIA = (
        "ligand_substruct",
        "ligand_substruct_pos",
        "ligand_max_atoms",
        "ligand_min_molweight",
        "ligand_max_molweight",
    )

    def has_rdkit_criteria(self) -> bool:
        """True if this leaf sets a criterion that needs an in-memory ligand scan."""
        return any(
            getattr(self, key) not in (None, [], "", 0) for key in self.RDKIT_CRITERIA
        )

    def has_criteria(self) -> bool:
        """True if this leaf sets any filter criterion at all."""
        return self.asdict() != _FILTER_DEFAULTS

    def rdkit_only(self) -> Union["Filter", None]:
        """A copy carrying only this leaf's RDKit criteria, or None if it has none.

        Used when flat SMARTS/property criteria are given alongside a nested expression:
        they apply to the whole expression, so they become their own leaf ANDed onto it.

        Returns:
            Filter | None: a leaf with just the RDKit criteria (plus ligand_operator,
            which controls how multiple SMARTS combine)
        """
        criteria = {
            key: getattr(self, key)
            for key in self.RDKIT_CRITERIA
            if getattr(self, key) not in (None, [], "", 0)
        }
        if not criteria:
            return None
        if self.ligand_operator:
            criteria["ligand_operator"] = self.ligand_operator
        return Filter(**criteria)

    @classmethod
    def get_filter_keys(cls, group: str) -> list:
        """Provide keys associated with each of the filter groups.
        Args:
            group (str): includese property filters, interaction filters, ligand filters, or all filters
        Returns:
            list of filter keywords associated with the specified group(s)
        """

        if group.lower() not in ["property", "interaction", "ligand", "all"]:
            raise OptionError(
                f'{group.lower()} is not a valid filter group. Please use "property", "interactions", "ligand", or "all'
            )

        filter_groups = {
            "property": [
                "eworst",
                "ebest",
                "leworst",
                "lebest",
                "score_percentile",
                "le_percentile",
            ],
            "interaction": [
                "vdw_interactions",
                "hb_interactions",
                "reactive_interactions",
            ],
            "ligand": [
                "ligand_name",
                "ligand_substruct",
                "ligand_substruct_pos",
                "ligand_max_atoms",
                "ligand_operator",
                "ligand_min_molweight",
                "ligand_max_molweight",
            ],
        }
        if group.lower() == "all":
            lst = []
            for i in filter_groups:
                lst.extend(filter_groups[i])
            return lst
        else:
            lst = filter_groups[group.lower()]
        return lst


# an all-unset leaf, for telling "no criteria" from "some criteria"
_FILTER_DEFAULTS = Filter().asdict()


class Filters:
    """A whole filter specification: a boolean tree of AND/OR nodes over :class:`Filter`
    leaves.

    A flat filter set is just a tree with one leaf, so the two ways of filtering share one
    type. Construct it either way::

        Filters(eworst=-9, hb_interactions=[("A:VAL:279:", True)])   # one leaf
        Filters.from_dict({"op": "or", "children": [{...}, {...}]})  # nested

    Criteria within a leaf are combined with AND; leaves are combined by their parent
    node's operator. Knowing whether a node is a boolean node or a leaf lives here, so
    callers traverse with ``leaves()`` / ``children`` instead of inspecting dicts.
    """

    OPERATORS = ("and", "or")

    def __init__(self, op: str = "and", children: list = None, **leaf_criteria):
        """
        Args:
            op (str): "and" or "or", how this node's children combine
            children (list, optional): child nodes, each a Filter or a Filters
            **leaf_criteria: any Filter field; given these, this becomes a single-leaf
                node holding one Filter built from them (the flat case)

        Raises:
            OptionError: on an unknown operator, or if both children and leaf criteria
                are given
        """
        op = str(op).strip().lower()
        if op not in self.OPERATORS:
            raise OptionError(
                f"Unrecognized boolean operator in filters: {op!r}. Use 'and' or 'or'."
            )
        if children and leaf_criteria:
            raise OptionError(
                "Filters takes either child nodes or flat filter criteria, not both."
            )
        self.op = op
        if leaf_criteria:
            self.children = [Filter(**leaf_criteria)]
        else:
            self.children = list(children) if children else []

    @classmethod
    def from_dict(cls, data) -> "Filters":
        """Build a filter tree from a plain dict (the JSON/GUI wire format).

        Accepts either a boolean node ``{"op": ..., "children": [...]}`` or a bare leaf
        dict of Filter fields — a bare leaf becomes a single-leaf tree, which is also how
        filter specifications saved by older versions are read back.

        Args:
            data (dict | Filters | Filter): the expression to parse

        Returns:
            Filters: the parsed tree

        Raises:
            OptionError: on an unknown operator or a malformed node
        """
        if isinstance(data, cls):
            return data
        if isinstance(data, Filter):
            return cls(children=[data])
        if not isinstance(data, dict):
            raise OptionError(
                f"Cannot build filters from {type(data).__name__}; expected a dict."
            )
        if "op" in data or "children" in data:
            if "children" not in data:
                raise OptionError("Filter node has an operator but no 'children'.")
            children = data["children"]
            if not isinstance(children, list):
                raise OptionError("Filter node 'children' must be a list.")
            return cls(
                op=data.get("op", "and"),
                children=[cls._child_from_dict(child) for child in children],
            )
        return cls(children=[Filter(**data)])

    @classmethod
    def _child_from_dict(cls, child):
        """A child is another boolean node, or a leaf dict of Filter fields."""
        if isinstance(child, (cls, Filter)):
            return child
        if not isinstance(child, dict):
            raise OptionError(
                f"Filter node children must be dicts, got {type(child).__name__}."
            )
        if "op" in child or "children" in child:
            return cls.from_dict(child)
        return Filter(**child)

    def to_dict(self) -> dict:
        """Serialize to plain nested dicts, ready for JSON.

        A single-leaf tree serializes to the bare leaf dict rather than a node wrapping
        one child — the two mean the same thing, and it keeps flat filter specifications
        byte-identical to what earlier versions wrote. ``from_dict`` round-trips both.

        Returns:
            dict: the leaf's fields if flat, else ``{"op": ..., "children": [...]}``
        """
        if self.is_flat():
            return self.leaf.asdict()
        return {
            "op": self.op,
            "children": [
                child.to_dict() if isinstance(child, Filters) else child.asdict()
                for child in self.children
            ],
        }

    def is_flat(self) -> bool:
        """True if this is a plain filter set — exactly one leaf and no nesting."""
        return len(self.children) == 1 and isinstance(self.children[0], Filter)

    @property
    def leaf(self) -> Filter:
        """The single Filter of a flat specification.

        Raises:
            OptionError: if this is a nested tree, where no single leaf applies
        """
        if not self.is_flat():
            raise OptionError(
                "These filters are a nested expression and have no single leaf; "
                "iterate leaves() instead."
            )
        return self.children[0]

    def leaves(self):
        """Yield every Filter in the tree, depth first."""
        for child in self.children:
            if isinstance(child, Filters):
                yield from child.leaves()
            else:
                yield child

    def is_empty(self) -> bool:
        """True if no leaf sets any criterion."""
        return not any(leaf.has_criteria() for leaf in self.leaves())

    def uses_percentile(self) -> bool:
        """True if any leaf filters on a score or ligand-efficiency percentile."""
        return any(
            leaf.score_percentile is not None or leaf.le_percentile is not None
            for leaf in self.leaves()
        )

    def rdkit_group_count(self) -> int:
        """How many top-level groups contain an RDKit criterion.

        Each one costs its own in-memory ligand scan, so more than one is worth warning
        about.
        """
        count = 0
        for child in self.children:
            leaves = child.leaves() if isinstance(child, Filters) else [child]
            if any(leaf.has_rdkit_criteria() for leaf in leaves):
                count += 1
        return count

    get_filter_keys = Filter.get_filter_keys

    def __eq__(self, other):
        if not isinstance(other, Filters):
            return NotImplemented
        return self.op == other.op and self.children == other.children

    def __repr__(self):
        if self.is_flat():
            return f"Filters({self.leaf!r})"
        return f"Filters(op={self.op!r}, children={self.children!r})"
