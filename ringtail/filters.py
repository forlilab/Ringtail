#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail Filter dataclass
#

from dataclasses import dataclass, field, asdict, astuple
import typing
from .exceptions import OptionError


@dataclass
class Filters:
    # property filters
    eworst: float = None
    ebest: float = None
    leworst: float = None
    lebest: float = None
    score_percentile: float = None
    le_percentile: float = None
    # interaction filters
    vdw_interactions: list[tuple] = field(default_factory=list)  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
    hb_interactions: list[tuple] = field(default_factory=list)  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
    reactive_interactions: list[tuple] = field(default_factory=list)  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
    interactions_count: list[tuple] = field(default_factory=list)  # e.g. [('hb_count', 5)]
    react_any: bool = None
    max_miss: int = 0
    # ligand filters
    ligand_name: list[str] = field(default_factory=list)  # e.g. ["lig1", "lig2"]
    ligand_substruct: list[str] = field(default_factory=list)  # e.g. ["ccc", "CN"]
    ligand_substruct_pos: list[str] = field(default_factory=list)  # e.g. ['"[Oh]C" 0 1.2 -5.5 10.0 15.5'] -> ["smart_string index_of_positioned_atom cutoff_distance x y z"]
    ligand_max_atoms: int = None
    ligand_operator: str = "OR"  # Choose AND or OR
    filter_ligands_flag: bool = False

    def __post_init__(self):

        # check percentiles between 0 and 100
        if self.score_percentile is not None and (self.score_percentile < 0 or self.score_percentile > 100):
            raise OptionError(f"Given score_percentile {self.score_percentile} not allowed. Should be within percentile range of 0-100.")
        if self.le_percentile is not None and (self.le_percentile < 0 or self.le_percentile > 100):
            raise OptionError(f"Given score_percentile {self.le_percentile} not allowed. Should be within percentile range of 0-100.")

        # check ligand operator
        if self.ligand_operator not in ["OR", "AND"]:
            raise OptionError(f"Given ligand_operator {self.ligand_operator} not allowed. Must be 'OR' or 'AND'.")

    @classmethod
    def get_defaults(cls):
        return cls().__dict__

    @classmethod
    def get_default_types(cls):
        return typing.get_type_hints(cls.__init__)

    def to_dict(self):
        return asdict(self)

    def to_list(self):
        return list(astuple(self))
