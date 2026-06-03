#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail output manager
#

from collections.abc import Iterable
from .exceptions import OutputError
from .ringtailoptions import Filters
import os
import json
from rdkit import Chem
from rdkit.Chem import SDWriter
from meeko import RDKitMolCreate


def _mol_with_one_conformer(mol: Chem.Mol, conf_id: int) -> Chem.Mol:
    """Return a copy of mol retaining only the specified conformer."""
    copy = Chem.RWMol(mol)
    for cid in [c.GetId() for c in mol.GetConformers() if c.GetId() != conf_id]:
        copy.RemoveConformer(cid)
    return copy.GetMol()


class OutputManager:
    """Class for creating outputs, can be a context manager to handle log files

    Attributes:
        output_log (str): name for log file
        _log_open (bool): if log file is open or not
    """

    def __init__(self, output_log: str = "output_log.txt", append_to_file: bool = True):
        self.output_log = output_log
        self._log_file = None
        self._log_open = False
        self.append = append_to_file

    def __enter__(self):
        """Opening outputmanager as a context manager"""
        self.open_logfile()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Closing outputmanager as a context manager"""
        self.close_logfile()
        return False

    # -#-#- Log file methods -#-#-#
    def open_logfile(self, write_filters_header=True):
        """
        Opens log file and creates it if needed

        Args:
            write_filters_header (bool): only used because one method does not take the same headers

        Raises:
            OutputError
        """
        if self.append:
            open_mode = "a"
        else:
            open_mode = "w"
        try:
            self._log_file = open(self.output_log, open_mode)
            self._log_open = True
            if write_filters_header:
                self._log_file.write("Filters:\n")
                self._log_file.write("***************\n")
        except Exception as e:
            raise OutputError("Error while creating log file") from e

    def close_logfile(self):
        """Closes the log file."""
        if self._log_open:
            self._log_file.close()
            self._log_file = None
            self._log_open = False

    def write_filter_results_in_log(self, lines: Iterable[tuple]):
        """Writes lines from results iterable into log file

        Args:
            lines (iterable): Iterable with tuples of data for
                writing into log

        Raises:
            OutputError
        """
        try:
            formatted_lines = []
            for line in lines:
                formatted_tuple = tuple(
                    round(item, 2) if isinstance(item, float) else item for item in line
                )
                formatted_lines.append(
                    ", ".join(
                        str(item) if not isinstance(item, str) else item
                        for item in formatted_tuple
                    )
                )
            self._log_file.write("\n".join(formatted_lines) + "\n***************\n\n")
        except Exception as e:
            raise OutputError("Error occurred during log writing") from e

    def _write_log_line(self, line: str):
        """write a single row to the log file

        Args:
            line (str): Line to write to log

        Raises:
            OutputError
        """
        try:
            self._log_file.write(line)
            self._log_file.write("\n")
        except Exception as e:
            raise OutputError(f"Error writing line {line} to log") from e

    def log_num_passing_ligands(self, number_passing_ligands: int):
        """
        Write the number of ligands which pass given filter to log file

        Args:
            number_passing_ligands (int): number of ligands that passed filter

        Raises:
            OutputError
        """
        try:
            self._log_file.write(
                f"\nNumber passing ligands: {number_passing_ligands} \n---------------\n"
            )
        except Exception as e:
            raise OutputError("Error writing number of passing ligands in log") from e

    def write_bookmarkname_in_log(self, bookmark_name):
        """Write the name of the result bookmark into log

        Args:
            bookmark_name (str): name of current results' bookmark in db

        Raises:
            OutputError
        """
        try:
            self._log_file.write(
                f"\nResult bookmark name: {bookmark_name}\n***************\n"
            )
        except Exception as e:
            raise OutputError("Error writing bookmark name to log") from e

    def write_filtervalues_in_log(
        self,
        filters_dict: dict,
        included_interactions: list,
        bookmark_name: str,
        additional_info="",
    ):
        """Takes dictionary of filters, formats as string and writes to log file

        Args:
            filters_dict (dict): dictionary with filtering options
            included_interactions (list): types of interactions to include in the filtering
            additional_info (str): any additional information to write to top of log file

        Raises:
            OutputError
        """
        try:
            buff = [additional_info, "##### PROPERTIES"]
            filters = filters_dict.copy()
            for k in Filters.get_filter_keys("property"):
                v = filters.pop(k, None)
                if v is not None:
                    v = "%2.3f" % v
                else:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))
            buff.append("#### LIGAND FILTERS")
            for k in Filters.get_filter_keys("ligand"):
                v = filters.pop(k, None)
                if v is not None:
                    if isinstance(v, list):
                        v = ", ".join([str(f) for f in v if f != ""])
                else:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))
            buff.append("#### INTERACTIONS")
            for _type in Filters.get_filter_keys("interaction"):
                info = filters.pop(_type, None)
                kept_interactions = []
                if len(info or []) == 0:
                    buff.append("#  % 7s :  [ none ]" % (_type))
                    continue
                for interact in info:
                    if _type + "-" + interact[0] not in included_interactions:
                        continue
                    else:
                        kept_interactions.append(interact)
                res_str = ", ".join(
                    ["(%s)%s" % ("~" if x[1] else "", x[0]) for x in kept_interactions]
                )
                label_str = "#  % 7s : %s" % (_type, res_str)
                buff.append(label_str)

            buff.append("#### OTHER FILTERS")
            for k, v in filters.items():
                if v is None:
                    v = " [ none ]"
                buff.append("#  % 7s : %s" % (k, v))

            for line in buff:
                self._write_log_line(line)

            self.write_bookmarkname_in_log(bookmark_name)

        except Exception as e:
            raise OutputError("Error occurred while writing filters to log") from e

    def write_maxmiss_union_header(self):
        """
        Properly formats header for the log file if using max_miss and enumerate_interaction_combs
        """
        self._log_file.write("\n---------------\n")
        self._log_file.write("Max Miss Union:\n")

    def write_find_similar_header(self, query_ligname, cluster_name):
        """
        Properly formats header for the log file find_similar_ligands
        """
        if not self._log_open:
            self.open_logfile(write_filters_header=False)
        self._log_file.write("\n---------------\n")
        self._log_file.write(
            f"Found ligands similar to {query_ligname} in clustering {cluster_name}:\n"
        )

    # -#-#- Non-logfile methods -#-#-#

    @staticmethod
    def _serialize_properties(properties: dict) -> dict:
        return {
            k: json.dumps(v) if isinstance(v, list) else v if isinstance(v, str) else str(v)
            for k, v in properties.items()
        }

    def write_out_mols(self, filename: str, mol_entries, export_sdf_directory: str):
        """Batch-write multiple mols to one SDF file with a single file-open.

        Args:
            filename (str): SDF filename
            mol_entries: iterable of (mol, flexres_mols, properties)
            export_sdf_directory (str): directory to write into
        """
        filepath = os.path.join(export_sdf_directory, filename)
        try:
            with open(filepath, "w") as sdf_file:
                w = SDWriter(sdf_file)
                for mol, flexres_mols, properties in mol_entries:
                    combined = RDKitMolCreate.combine_rdkit_mols([mol] + flexres_mols)
                    for k, v in self._serialize_properties(properties).items():
                        combined.SetProp(k, v)
                    for conf in combined.GetConformers():
                        w.write(combined, conf.GetId())
                w.close()
        except Exception as e:
            raise OutputError("Error occurred while batch-writing SDF") from e

    def write_out_mol(
        self, filename, mol, flexres_per_conf, properties, export_sdf_directory: str
    ):
        """Writes out given mol as sdf. Will create the specified sdf folder in
        current working directory if needed.

        Args:
            filename (str): name of SDF file that will be written to
            mol (Chem.Mol): RDKit mol with 1+ conformers
            flexres_per_conf (list[list[Chem.Mol]]): per-conformer flexres mol lists,
                or [] if no flexres
            properties (dict): properties to embed in each SDF record

        Raises:
            OutputError
        """
        filepath = os.path.join(export_sdf_directory, filename)
        try:
            str_props = self._serialize_properties(properties)
            with open(filepath, "a") as sdf_file:
                w = SDWriter(sdf_file)
                if not flexres_per_conf:
                    for k, v in str_props.items():
                        mol.SetProp(k, v)
                    for conf in mol.GetConformers():
                        w.write(mol, conf.GetId())
                else:
                    for conf_idx, conf in enumerate(mol.GetConformers()):
                        single = _mol_with_one_conformer(mol, conf.GetId())
                        flexres = (
                            flexres_per_conf[conf_idx]
                            if conf_idx < len(flexres_per_conf)
                            else []
                        )
                        combined = RDKitMolCreate.combine_rdkit_mols([single] + flexres)
                        for k, v in str_props.items():
                            combined.SetProp(k, v)
                        w.write(combined, combined.GetConformer().GetId())
                w.close()
        except Exception as e:
            raise OutputError("Error occurred while writing SDF from RDKit Mol") from e

    def write_receptor_pdbqt(self, recname: str, receptor_str: str):
        """
        Writes a pdbqt file from receptor "blob"

        Args:
            recname (str): name of receptor to use in output filename
            receptor_compbytes (blob): receptor blob
        """

        try:
            if not recname.endswith(".pdbqt"):
                recname = recname + ".pdbqt"
            with open(recname, "w") as f:
                f.write(receptor_str)
        except Exception as e:
            raise OutputError(f"Error writing receptor pdbqt to {recname}") from e
