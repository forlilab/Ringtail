#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail Filter dataclass
#

from dataclasses import dataclass, field, fields, asdict, astuple
import typing
from .exceptions import OptionError
from .logmanager import logger

class TypeSafe:
    """
    Class that handles safe typesetting of values in the other option classes. 
    It creates internal attribtues that are get/set according to the type that 
    is specified. Raises error if wrong type is attempted. 
    """
    def __init__(self, attr, attrtype):                        
        self.attrpublic = attr
        self.attrprivate = "_" + attr
        self.type = attrtype

    def __get__(self, obj, objtype=None):  
        value = getattr(obj, self.attrprivate)
        return value
    
    def __set__(self, obj, value):
        if type(value) in (self.type, type(None)):
            setattr(obj, self.attrprivate, value) 
            logger.info(f'{self.attrpublic} was set to {value}.')
        else:
            raise OptionError(f'{self.attrpublic} can only be of type {self.type}, but was attempted set as {type(value)} which is invalid.')

class Filters:
    """
    Object that holds all optional filters.
    
    Args:
        eworst (float): 
        ebest (float): 
        leworst (float): 
        score_percentile (float):
        le_percentile (float):
        vdw_interactions (list[tuple]): e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
        hb_interactions (list[tuple]): e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
        reactive_interactions (list[tuple]): e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
        interactions_count (list[tuple]): e.g. [('hb_count', 5)]
        react_any (bool): 
        max_miss (int): 
        ligand_name (list[str]): e.g. ["lig1", "lig2"]
        ligand_substruct (list[str]): e.g. ["ccc", "CN"]
        ligand_substruct_pos (list[str]): e.g. ['"[Oh]C" 0 1.2 -5.5 10.0 15.5'] -> ["smart_string index_of_positioned_atom cutoff_distance x y z"]
        ligand_max_atoms (int): 
        ligand_operator (str): AND or OR
        """
    
    eworst = TypeSafe("eworst", float)
    ebest = TypeSafe("ebest", float)
    leworst = TypeSafe("leworst", float)
    lebest = TypeSafe("lebest", float)
    score_percentile = TypeSafe("score_percentile", float)
    le_percentile = TypeSafe("le_percentile", float)
    vdw_interactions = TypeSafe("wdv_interactions", list) #of tuple
    hb_interactions = TypeSafe("hb_interactions", list) #of tuple
    reactive_interactions = TypeSafe("reactive_interactions", list) #of tuple
    interactions_count = TypeSafe("interactions_count", list) #of tuple
    react_any = TypeSafe("react_any", bool)
    max_miss = TypeSafe("max_miss", int)
    ligand_name = TypeSafe("ligand_name", list) #of str
    ligand_substruct = TypeSafe("ligand_substruct", list) #of str
    ligand_substruct_pos = TypeSafe("ligand_substruct_pos", list) #of str
    ligand_max_atoms = TypeSafe("ligand_max_atoms", int)
    ligand_operator = TypeSafe("ligand_operator", str)

    def __init__(self,
                 eworst = None,
                 ebest = None,
                 leworst = None,
                 lebest = None,
                 score_percentile = None,
                 le_percentile = None,
                 vdw_interactions = [],
                 hb_interactions = [],
                 reactive_interactions = [],
                 interactions_count = [],
                 react_any = None,
                 max_miss = 0,
                 ligand_name = [],
                 ligand_substruct = [],
                 ligand_substruct_pos = [],
                 ligand_max_atoms = None,
                 ligand_operator = "OR"):

        self.eworst = eworst
        self.ebest = ebest
        self.leworst = leworst
        self.lebest = lebest
        self.score_percentile = score_percentile
        self.le_percentile = le_percentile
        self.vdw_interactions = vdw_interactions
        self.hb_interactions = hb_interactions
        self.reactive_interactions = reactive_interactions
        self.interactions_count = interactions_count
        self.react_any = react_any
        self.max_miss = max_miss
        self.ligand_name = ligand_name
        self.ligand_substruct = ligand_substruct
        self.ligand_substruct_pos = ligand_substruct_pos
        self.ligand_max_atoms = ligand_max_atoms
        self.ligand_operator = ligand_operator
    
    def __setattr__(self, name, value):
        """ Overloaded setattr method so value compatibility check can be ran"""
        super().__setattr__(name, value)
        if hasattr(self, "ligand_operator"):
            # ligand_operator is value to be set during init, ensures all values present before first check
            self._compatibility_checks()
        
    def _compatibility_checks(self):
        """
        Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed.
        
        #TODO:
            - some of these should probably be warnings, and maintain the object but giving you a chance to correct the bad filter
        """
        
        if (self.eworst is not None and self.score_percentile is not None):
            logger.warning("Cannot use --eworst cutoff with --score_percentile. Overiding score_percentile with eworst.")
            self.score_percentile = None
        
        if (self.leworst is not None and self.le_percentile is not None):
            logger.warning("Cannot use --leworst cutoff with --le_percentile. Overiding le_percentile with leworst.")
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
        """
        Makes a dict of filter kv based on given group. 
        Args:
            group (str): either "property", "interactions", or "ligand"
        """
        if group.lower() not in ["property", "interaction", "ligand"]:
            raise OptionError(f'{group.lower()} is not a valid filter group. Please use "property", "interactions", or "ligand"')

        filter_groups = {
        "property": [
            "eworst",
            "ebest",
            "leworst",
            "lebest",
            "score_percentile",
            "le_percentile"
            ],
        "interaction": [
            "vdw_interactions",
            "hb_interactions",
            "reactive_interactions"
        ],
        "ligand": [
            "ligand_name",
            "ligand_substruct",
            "ligand_substruct_pos",
            "ligand_max_atoms",
            "ligand_operator"
        ]
    }
        
        filters = {}
        for i in filter_groups[group.lower()]:
            filters[i] = getattr(self, i)

        return filters

    def todict(self) -> dict:
        returndict = {}
        attributes = [attr for attr in vars(self) 
                if (not attr.startswith('__')
                )]
        for attribute in attributes:
            returndict[attribute.strip('_')] = getattr(self, attribute)
        return returndict
# @dataclass
# class Filters:
#     # property filters
#     eworst: float = field(default=None, metadata={"filter_type": "property"})
#     ebest: float = field(default=None, metadata={"filter_type": "property"})
#     leworst: float = field(default=None, metadata={"filter_type": "property"})
#     lebest: float = field(default=None, metadata={"filter_type": "property"})
#     score_percentile: float = field(default=None, metadata={"filter_type": "property"})
#     le_percentile: float = field(default=None, metadata={"filter_type": "property"})
#     # interaction filters
#     vdw_interactions: list[tuple] = field(default_factory=list, metadata={"filter_type": "interaction"})  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
#     hb_interactions: list[tuple] = field(default_factory=list, metadata={"filter_type": "interaction"})  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]
#     reactive_interactions: list[tuple] = field(default_factory=list, metadata={"filter_type": "interaction"})  # e.g. [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]

#     interactions_count: list[tuple] = field(default_factory=list, metadata={"filter_type": "interactions_count"})  # e.g. [('hb_count', 5)]
#     react_any: bool = field(default=None, metadata={"filter_type": "react_any"})
#     max_miss: int = field(default=0, metadata={"filter_type": "max_miss"})
#     # ligand filters
#     ligand_name: list[str] = field(default_factory=list, metadata={"filter_type": "ligand"})  # e.g. ["lig1", "lig2"]
#     ligand_substruct: list[str] = field(default_factory=list, metadata={"filter_type": "ligand"})  # e.g. ["ccc", "CN"]
#     ligand_substruct_pos: list[str] = field(default_factory=list, metadata={"filter_type": "ligand"})  # e.g. ['"[Oh]C" 0 1.2 -5.5 10.0 15.5'] -> ["smart_string index_of_positioned_atom cutoff_distance x y z"]
#     ligand_max_atoms: int = field(default=None, metadata={"filter_type": "ligand"})
#     ligand_operator: str = field(default="OR", metadata={"filter_type": "ligand"})  # Choose AND or OR
        
#     def __setattr__(self, name, value):
#         super().__setattr__(name, value)
#         self._compatibility_checks()
        
#     def _compatibility_checks(self):
#         '''It is probably not super efficient to run comp checks every time a value is set, but in this case attributes are not that many 
#         and is only set a limited number of times per ringtail core. I don't think it is a performance issue, but maybe a performance curiosity.'''

#         if (self.eworst is not None and self.score_percentile is not None):
#             logger.warning("Cannot use --eworst cutoff with --score_percentile. Overiding score_percentile with eworst.")
#             self.score_percentile = None
        
#         if (self.leworst is not None and self.le_percentile is not None):
#             logger.warning("Cannot use --leworst cutoff with --le_percentile. Overiding le_percentile with leworst.")
#             self.le_percentile = None  

#         if self.score_percentile is not None and (self.score_percentile < 0 or self.score_percentile > 100):
#             raise OptionError(f"Given score_percentile {self.score_percentile} not allowed. Should be within percentile range of 0-100.")
        
#         if self.le_percentile is not None and (self.le_percentile < 0 or self.le_percentile > 100):
#             raise OptionError(f"Given score_percentile {self.le_percentile} not allowed. Should be within percentile range of 0-100.")

#         if self.ligand_operator not in ["OR", "AND"]:
#             raise OptionError(f"Given ligand_operator {self.ligand_operator} not allowed. Must be 'OR' or 'AND'.")
        
#         if self.max_miss < 0:
#             raise OptionError("--max_miss must be greater than or equal to 0")

#     def filters_in_group(self, group: str) -> dict:
#         ''' Returns filters as a dict for a given group'''
#         fg = {}
#         for f in fields(self):
#             if f.metadata["filter_type"] == group:
#                 fg[f.name] = getattr(self, f.name)

#         return fg

#     @classmethod
#     def get_interaction_filter_keys(cls):
#         keys = []
#         for f in fields(cls):
#             if f.metadata["filter_type"] == "interaction":
#                 keys.append(f.name)
#         return keys

#     @classmethod
#     def get_property_filter_keys(cls):
#         keys = []
#         for f in fields(cls):
#             if f.metadata["filter_type"] == "property":
#                 keys.append(f.name)
#         return keys

#     @classmethod
#     def get_ligand_filter_keys(cls):
#         keys = []
#         for f in fields(cls):
#             if f.metadata["filter_type"] == "ligand":
#                 keys.append(f.name)
#         return keys

#     def todict(self):
#         return asdict(self)

#     def tolist(self):
#         return list(astuple(self))
