import os

from .exceptions import OptionError
from .logmanager import logger
import copy
import json

""" Ringtail options contains objects for holding all ringtail options, 
and ensures safe type enforcement."""

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
            logger.debug(f'{self.attrpublic} was set to {value}.')
        else:
            raise OptionError(f'{self.attrpublic} can only be of type {self.type}, but was attempted set as {type(value)} which is invalid.')

class NumberSafe(TypeSafe):
    def __init__(self, attr, attrtype):
        super().__init__(attr, attrtype)

    def __set__(self, obj, value):
        if type(value) in (float, int, type(None)):
            setattr(obj, self.attrprivate, value) 
            logger.debug(f'{self.attrpublic} was set to {value}.')
        else:
            raise OptionError(f'{self.attrpublic} can only be of type {self.type}, but was attempted set as {type(value)} which is invalid.')

class DirSafe(TypeSafe):
    def __init__(self, attr, attrtype):
        super().__init__(attr, attrtype)

    def __set__(self, obj, value: str):
        if type(value) in (str, type(None)):
            if value is not None and not value == "" and not value.endswith("/"): value += "/"
            setattr(obj, self.attrprivate, value) 
            logger.debug(f'{self.attrpublic} was set to {value}.')
        else:
            raise OptionError(f'{self.attrpublic} can only be of type {self.type}, but was attempted set as {type(value)} which is invalid.')

class RTOptions:
    """ Holds standard methods for the ringtail option child classes."""
    @classmethod
    def is_valid_path(self, path):
        if os.path.exists(path):
            return True
        else:
            return False
    
    #TODO method to get values types names/keywords etc
        
    @classmethod
    def get_all_options_types(self):
        """
        This needs to be a method that grabs attributes from each class, and figures out what type they are as a TypeSafe object
        Also grab their default values
        """
        pass

    def todict(self) -> dict:
        dict = {}
        attributes = [attr for attr in vars(self) 
                if (not attr.startswith('__')
                )]
        for attribute in attributes:
            dict[attribute.strip('_')] = getattr(self, attribute)
        returndict = copy.deepcopy(dict)
        
        return returndict

class GeneralOptions(RTOptions):
    """Creates a class with general options relevant to any Ringtail process, including:
    Args:
        docking_mode: specify AutoDock program used to generate results. Available options are "DLG" and "Vina". Vina mode will automatically change --pattern to *.pdbqt
        summary: prints summary information about stored data to STDOUT.
        verbose: Print results passing filtering criteria to STDOUT. NOTE: runtime may be slower option used
        debug: Print additional error information to STDOUT
        """
    options = {
        "docking_mode":{
            "type": str,
            "default": "dlg"
        },
        "summary":{
            "type": bool,
            "default": False
        },  
        "verbose":{
            "type": bool,
            "default": False
        },  
        "debug":{
            "type": bool,
            "default": False
        },  
    }
    
    docking_mode = TypeSafe("docking_mode", str)
    summary = TypeSafe("summary", bool)
    verbose = TypeSafe("verbose", bool)
    debug = TypeSafe("debug", bool)

    def __init__(self,
                 docking_mode="dlg",
                 summary=False,
                 verbose=False,
                 debug=False,):              
        self.docking_mode = docking_mode                        
        self.summary = summary                          
        self.verbose = verbose         
        self.debug = debug

class InputFiles(RTOptions):
    """ Class that handles sources of data to be written including ligand data paths and how 
    to traverse them, and options to store receptor.
    
    Args:
        file(list): ligand docking output file to save. Compressed (.gz) files allowed. Only 1 receptor allowed.
        file_path (list): directory(s) containing docking output files to save. Compressed (.gz) files allowed
        file_list (list): file(s) containing the list of docking output files to save; relative or absolute paths are allowed. Compressed (.gz) files allowed
        file_pattern (str): specify which pattern to use when searching for result files to process [only with "--file_path"]
        recursive (bool): enable recursive directory scan when --file_path is used
        receptor_file (str): Use with Vina mode. Give file for receptor PDBQT.
        save_receptor (bool): Saves receptor PDBQT to database. Receptor location must be specied with in --file, --file_path directory or --file_list file
        target (str): name of receptor
        """
    
    file=TypeSafe("file", list)
    file_path=TypeSafe("file_path", list)
    file_list=TypeSafe("file_list", list)
    file_pattern=TypeSafe("file_pattern", str)
    recursive=TypeSafe("recursive", bool)
    receptor_file=TypeSafe("receptor_file", str)
    save_receptor=TypeSafe("save_receptor", bool)
    target = None
    
    def __init__(self,
                 file=None, 
                 file_path=None, 
                 file_list=None, 
                 file_pattern='*.dlg*', 
                 recursive=None, 
                 receptor_file=None,
                 save_receptor=None,):
        if receptor_file is None:
            pass
        elif receptor_file is not None and self.is_valid_path(receptor_file):
            self.target = (os.path.basename(receptor_file).split(".")[0])
        else:
            raise OptionError("The receptor PDBQT file path is not valid. Please check location of receptor file and --receptor_file option")
        self.file=file 
        self.file_path=file_path 
        self.file_list=file_list
        self.file_pattern=file_pattern
        self.recursive=recursive
        self.receptor_file=receptor_file
        self.save_receptor=save_receptor
    
    #TODO on change need to add files
        
class ResultsProcessingOptions(RTOptions):
    """ Class that holds database write options that affects write time, such as how to 
    break up data files, number of computer processes to use, and and how many poses to store.
    
    Args:
        store_all_poses (bool): Store all poses from input files. Overrides --max_poses
        max_poses (int): n: Store top pose for top n clusters
        add_interactions (bool): Find interactions between ligand poses and receptor and save to database. Requires receptor PDBQT to be given with input files (all modes) and --receptor_file to be specified with Vina mode. SIGNIFICANTLY INCREASES DATBASE WRITE TIME.
        interaction_tolerance (float): Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired. Only compatible with ADGPU mode
        interaction_cutoffs (list): Use with --add_interactions, specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. These are the default cutoffs.
        max_proc (int): Maximum number of processes to create during parallel file parsing. Defaults to number of CPU processors.
        
        """

    store_all_poses = TypeSafe("store_all_poses", bool)
    max_poses = TypeSafe("max_poses", int)
    add_interactions = TypeSafe("add_interactions", bool)
    interaction_tolerance = NumberSafe("interaction_tolerance", float)
    interaction_cutoffs = TypeSafe("interaction_cutoffs", list)
    max_proc = TypeSafe("max_proc", int)

    def __init__(self,
                 store_all_poses = False,
                 max_poses = 3,
                 add_interactions = False,
                 interaction_tolerance = None,
                 interaction_cutoffs = [3.7, 4.0],
                 max_proc = None,):

        # Check interaction_tolerance compatibilities 
        if store_all_poses is True and interaction_tolerance is not None:
            logger.warning("Cannot use interaction_tolerance with store_all_poses. Removing interaction_tolerance.")
            self.interaction_tolerance = None
        else:
            self.interaction_tolerance = interaction_tolerance            
        
        self.store_all_poses = store_all_poses
        self.max_poses = max_poses
        self.add_interactions = add_interactions
        self.interaction_cutoffs = interaction_cutoffs
        self.max_proc = max_proc

    #TODO need onchange

class StorageOptions(RTOptions):
    """ Class that handles options for the storage (database) manager class, including
    conflict handling and result clustering and ordering.

    Args:
    filter_bookmark (str): Perform filtering over specified bookmark. (in output group in CLI)
    append_results (bool): Add new results to an existing database, specified by database choice in ringtail initialization or --input_db in cli
    duplicate_handling (str, options): specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.
    overwrite_logfile (bool): by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database
    order_results_by (str): Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION."
            "available fields are:  "
            '"e" (docking_score), '
            '"le" (ligand efficiency), '
            '"delta" (delta energy from best pose), '
            '"ref_rmsd" (RMSD to reference pose), '
            '"e_inter" (intermolecular energy), '
            '"e_vdw" (van der waals energy), '
            '"e_elec" (electrostatic energy), '
            '"e_intra" (intermolecular energy), '
            '"n_interact" (number of interactions), '
            '"rank" (rank of ligand pose), '
            '"run" (run number for ligand pose), '
            '"hb" (hydrogen bonds); '
    outfields (str): defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are:  '
            '"Ligand_name" (Ligand name), '
            '"e" (docking_score), '
            '"le" (ligand efficiency), '
            '"delta" (delta energy from best pose), '
            '"ref_rmsd" (RMSD to reference pose), '
            '"e_inter" (intermolecular energy), '
            '"e_vdw" (van der waals energy), '
            '"e_elec" (electrostatic energy), '
            '"e_intra" (intermolecular energy), '
            '"n_interact" (number of interactions), '
            '"ligand_smile" , '
            '"rank" (rank of ligand pose), '
            '"run" (run number for ligand pose), '
            '"hb" (hydrogen bonds), '
            '"receptor" (receptor name); '
            "Fields are printed in the order in which they are provided. Ligand name will always be returned and will be added in first position if not specified.
    output_all_poses (bool): By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.
    mfpt_cluster (float): Cluster filered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.
    interaction_cluster (float): Cluster filered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.
    results_view_name (str): name for resulting book mark file. Default value is "passing_results"

    """

    filter_bookmark = TypeSafe("filter_bookmark", str)
    append_results = TypeSafe("append_results", bool)
    duplicate_handling = TypeSafe("duplicate_handling", str) 
    overwrite = TypeSafe("overwrite", bool)
    order_results = TypeSafe("order_results", str)
    outfields = TypeSafe("outfields", str)
    output_all_poses = TypeSafe("output_all_poses", bool)
    mfpt_cluster = NumberSafe("mfpt_cluster", float)
    interaction_cluster = NumberSafe("interaction_cluster", float)
    results_view_name = TypeSafe("results_view_name", str)

    def __init__(self,
                 filter_bookmark = None,
                 append_results = False,
                 duplicate_handling = None,
                 overwrite_log_file = None,
                 order_results_by = None,
                 outfields = "Ligand_name,e",
                 output_all_poses = None,
                 mfpt_cluster = None,
                 interaction_cluster = None,
                 results_view_name = "passing_results",):
        self.filter_bookmark = filter_bookmark
        self.append_results = append_results
        self.order_results = order_results_by
        self.overwrite = overwrite_log_file
        self.outfields = outfields
        self.output_all_poses = output_all_poses
        self.mfpt_cluster = mfpt_cluster
        self.interaction_cluster = interaction_cluster
        self.results_view_name = results_view_name

        if duplicate_handling is not None:
            self.duplicate_handling = duplicate_handling.upper()
            if self.duplicate_handling not in {"IGNORE", "REPLACE"}:
                logger.warning(
                    f"--duplicate_handling option {duplicate_handling} not allowed. Reverting to default behavior."
                )
                self.duplicate_handling = None
        else: self.duplicate_handling = duplicate_handling

        if self.order_results is not None and self.order_results not in self.order_options:
            raise OptionError(
                "Requested ording option that is not available. Please see --help for available options."
            )
        
    order_options = {
            "e",
            "le",
            "delta",
            "ref_rmsd",
            "e_inter",
            "e_vdw",
            "e_elec",
            "e_intra",
            "n_interact",
            "rank",
            "run",
            "hb",
        }
    #TODO on change need to add files

class ReadOptions(RTOptions):
    """ Class that holds options related to reading from the database, including format for
    result export and alternate ways of displaying the data (plotting),

    Args:
        filtering (bool): implicit argument set if there are any optional filters present
        plot (bool): Makes scatterplot of LE vs Best Energy, saves as scatter.png.
        find_similar_ligands (str): Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering.
        export_bookmark_csv (str): Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables
        export_bookmark_db (bool): Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db
        export_query_csv (str): Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]
        export_receptor (bool): Export stored receptor pdbqt. Will write to current directory.
        data_from_bookmark (bool): Write log of --outfields data for bookmark specified by --bookmark_name. Must use without any filters.
        pymol (bool): Lauch PyMOL session and plot of ligand efficiency vs docking score for molecules in bookmark specified with --bookmark_name. Will display molecule in PyMOL when clicked on plot. Will also open receptor if given.
        enumerate_interaction_combs (bool): #TODO
        log_file (str): by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file
        export_sdf_path (str): specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.
    """

    filtering = TypeSafe("filtering", bool)
    plot = TypeSafe("plot", bool)
    find_similar_ligands = TypeSafe("find_similar_ligand", str)
    export_bookmark_csv = TypeSafe("export_bookmark_csv", str)
    export_bookmark_db = TypeSafe("export_bookmark_db", bool)
    export_query_csv = TypeSafe("export_query_csv", str)
    export_receptor = TypeSafe("export_receptor", bool)
    data_from_bookmark = TypeSafe("data_from_bookmark", bool)
    pymol = TypeSafe("pymol", bool)
    enumerate_interaction_combs = TypeSafe("enumerate_interaction_combs", bool)
    log_file = TypeSafe("log_file", str)
    export_sdf_path = DirSafe("export_sdf_path", str)
    def __init__(self,
                 filtering = None,
                 plot = None,
                 find_similar_ligands = None,
                 export_bookmark_csv =  None,
                 export_bookmark_db =  None,
                 export_query_csv =  None,
                 export_receptor =  None,
                 data_from_bookmark =  None,
                 pymol =  None,
                 enumerate_interaction_combs =  None,
                 log_file="output_log.txt",
                 export_sdf_path=""):
        self.filtering = filtering
        self.plot = plot
        self.find_similar_ligands = find_similar_ligands
        self.export_bookmark_csv = export_bookmark_csv
        self.export_bookmark_db = export_bookmark_db
        self.export_query_csv = export_query_csv
        self.export_receptor = export_receptor
        self.data_from_bookmark = data_from_bookmark
        self.pymol = pymol
        self.enumerate_interaction_combs = enumerate_interaction_combs
        self.log_file = log_file
        self.export_sdf_path = export_sdf_path
    
class Filters(RTOptions):
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
    
    eworst = NumberSafe("eworst", float)
    ebest = NumberSafe("ebest", float)
    leworst = NumberSafe("leworst", float)
    lebest = NumberSafe("lebest", float)
    score_percentile = NumberSafe("score_percentile", float)
    le_percentile = NumberSafe("le_percentile", float)
    vdw_interactions = TypeSafe("vdw_interactions", list) #of tuple
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
    
    def __getattr__(self, obj, objtype=None):  
        privname = "_" + obj.__name__
        value = obj.__get__(obj, privname)
        return value
        
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
    
    @classmethod
    def get_filter_keys(self, group) -> list:
        if group.lower() not in ["property", "interaction", "ligand", "all"]:
            raise OptionError(f'{group.lower()} is not a valid filter group. Please use "property", "interactions", "ligand", or "all')

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
        ]}
        if group.lower() == "all":
            list = []
            for i in filter_groups:
                list.extend(filter_groups[i])
            return list
        else:
            list = filter_groups[group.lower()]
        return list
    