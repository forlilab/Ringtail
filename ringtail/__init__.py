#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail
#
import warnings

warnings.filterwarnings("ignore", "pkg_resources is deprecated", UserWarning, "prody")


from .ringtailcore import (
    RingtailCore,
    storage_types,
    get_valid_storageclass,
)
from .cloptionparser import CLOptionParser
from .logutils import LOGGER, get_logger, setup_logging
from .ringtailoptions import (
    RingtailDefaults,
    Filters,
    validate_file_pattern,
    ringtail_defaults,
)
from .storagemanager import StorageManager
from .storagemanager_sqlite import StorageManagerSQLite
from .storagemanager_duckdb import StorageManagerDuckDB
from .querybuilder import QueryBuilder
from .mpmanager import MPManager
from .mpreaderwriter import DockingFileReader, Writer
from .parsers import (
    VinaMoleculeSupplier,
    ADGPUMoleculeSupplier,
    SDFMoleculeSupplier,
)
from .receptormanager import ReceptorManager
from .resultsmanager import ResultsManager
from .outputmanager import OutputManager
from .interactions import InteractionFinder, find_interactions
from .clustermanager import (
    MorganFingerprintCluster,
    InteractionBitvectorCluster,
    top_score_per_cluster,
    butina_cluster_fingerprints,
)
from .exceptions import OptionError, RingtailError
from ._version import __version__

__all__ = [
    "RingtailCore",
    "__version__",
    "RingtailDefaults",
    "validate_file_pattern",
    "storage_types",
    "get_valid_storageclass",
    "Filters",
    "LOGGER",
    "get_logger",
    "setup_logging",
    "CLOptionParser",
    "StorageManager",
    "StorageManagerSQLite",
    "StorageManagerDuckDB",
    "QueryBuilder",
    "MPManager",
    "DockingFileReader",
    "Writer",
    "ADGPUMoleculeSupplier",
    "VinaMoleculeSupplier",
    "SDFMoleculeSupplier",
    "find_interactions",
    "ReceptorManager",
    "ResultsManager",
    "OutputManager",
    "InteractionFinder",
    "MorganFingerprintCluster",
    "InteractionBitvectorCluster",
    "top_score_per_cluster",
    "butina_cluster_fingerprints",
    "OptionError",
    "RingtailError",
    "ringtail_defaults",
]
