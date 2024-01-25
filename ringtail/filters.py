#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail Filter dataclass
#

from dataclasses import dataclass, field, fields, asdict, astuple
import typing
from .exceptions import OptionError
import logging


@dataclass
class Filters:
    # property filters
    eworst: float = field(default=None, metadata={"filter_type": "property"})
    ebest: float = field(default=None, metadata={"filter_type": "property"})
    leworst: float = field(default=None, metadata={"filter_type": "property"})
    lebest: float = field(default=None, metadata={"filter_type": "property"})
    score_percentile: float = field(default=None, metadata={"filter_type": "property"})
    le_percentile: float = field(default=None, metadata={"filter_type": "property"})
    # interaction filters
    vdw_interactions: list[tuple] = field(default_factory=list, metadata={"filter_type": "interaction"})  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
    hb_interactions: list[tuple] = field(default_factory=list, metadata={"filter_type": "interaction"})  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
    reactive_interactions: list[tuple] = field(default_factory=list, metadata={"filter_type": "interaction"})  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]

    interactions_count: list[tuple] = field(default_factory=list, metadata={"filter_type": "interactions_count"})  # e.g. [('hb_count', 5)]
    react_any: bool = field(default=None, metadata={"filter_type": "react_any"})
    max_miss: int = field(default=0, metadata={"filter_type": "max_miss"})
    # ligand filters
    ligand_name: list[str] = field(default_factory=list, metadata={"filter_type": "ligand"})  # e.g. ["lig1", "lig2"]
    ligand_substruct: list[str] = field(default_factory=list, metadata={"filter_type": "ligand"})  # e.g. ["ccc", "CN"]
    ligand_substruct_pos: list[str] = field(default_factory=list, metadata={"filter_type": "ligand"})  # e.g. ['"[Oh]C" 0 1.2 -5.5 10.0 15.5'] -> ["smart_string index_of_positioned_atom cutoff_distance x y z"]
    ligand_max_atoms: int = field(default=None, metadata={"filter_type": "ligand"})
    ligand_operator: str = field(default="OR", metadata={"filter_type": "ligand"})  # Choose AND or OR
        
    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        self._compatibility_checks()
        
    def _compatibility_checks(self):
        '''It is probably not super efficient to run comp checks every time a value is set, but in this case attributes are not that many 
        and is only set a limited number of times per ringtail core. I don't think it is a performance issue, but maybe a performance curiosity.'''

        if (self.eworst is not None and self.score_percentile is not None):
            logging.warning("Cannot use --eworst cutoff with --score_percentile. Overiding score_percentile with eworst.")
            self.score_percentile = None
        
        if (self.leworst is not None and self.le_percentile is not None):
            logging.warning("Cannot use --leworst cutoff with --le_percentile. Overiding le_percentile with leworst.")
            self.le_percentile = None  

        if self.score_percentile is not None and (self.score_percentile < 0 or self.score_percentile > 100):
            raise OptionError(f"Given score_percentile {self.score_percentile} not allowed. Should be within percentile range of 0-100.")
        
        if self.le_percentile is not None and (self.le_percentile < 0 or self.le_percentile > 100):
            raise OptionError(f"Given score_percentile {self.le_percentile} not allowed. Should be within percentile range of 0-100.")

        if self.ligand_operator not in ["OR", "AND"]:
            raise OptionError(f"Given ligand_operator {self.ligand_operator} not allowed. Must be 'OR' or 'AND'.")
        
        if self.max_miss < 0:
            raise OptionError("--max_miss must be greater than or equal to 0")

    def filters_in_group(self, group: str) -> dict:
        ''' Returns filters as a dict for a given group'''
        fg = {}
        for f in fields(self):
            if f.metadata["filter_type"] == group:
                fg[f.name] = getattr(self, f.name)

        return fg
        
    @classmethod
    def get_defaults(cls):
        return cls().__dict__

    @classmethod
    def get_default_types(cls):
        return typing.get_type_hints(cls.__init__)

    @classmethod
    def get_interaction_filter_keys(cls):
        keys = []
        for f in fields(cls):
            if f.metadata["filter_type"] == "interaction":
                keys.append(f.name)
        return keys

    @classmethod
    def get_property_filter_keys(cls):
        keys = []
        for f in fields(cls):
            if f.metadata["filter_type"] == "property":
                keys.append(f.name)
        return keys

    @classmethod
    def get_ligand_filter_keys(cls):
        keys = []
        for f in fields(cls):
            if f.metadata["filter_type"] == "ligand":
                keys.append(f.name)
        return keys

    def to_dict(self):
        return asdict(self)

    def to_list(self):
        return list(astuple(self))
