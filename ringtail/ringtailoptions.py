import os
from .exceptions import OptionError
from .logmanager import logger
import copy

""" Ringtail options contains objects for holding all ringtail options, 
and ensures safe type enforcement."""

class TypeSafe:
    """
    Class that handles safe typesetting of values of a specified built-in type. 
    Any TypeSafe object has a setter that checks the type, which is stored as an internal attribute 
    of the object. 
    It is the hope to extend this to work with custom types, such as "percentage" (float with a max and min value),
    and direcotry (string that must end with '/'). 
    Raises OptionError if wrong type is attempted.     
    """
    def __init__(self, default, type, object_name):
        self.object_name = object_name
        self.type = type
        self.default = default
        self.value = self.default
        
    
    def __setattr__(self, name, value):
        if name == "value":
            if type(value) == self.type:
                self.__dict__["value"] = value
                logger.debug(f"{self.object_name} updated to {self.value}.")
            elif self.type == float and type(value) in [float, int]: 
                self.__dict__["value"] = float(value)
                logger.debug(f"{self.object_name} updated to {self.value}.")
            else:
                self.__dict__["value"] = self.default
                logger.debug(f"{self.object_name} reset to default {self.default}.")
        else:
            self.__dict__[name] = value

class RTOptions:
    """ Holds standard methods for the ringtail option child classes.
    Options can be added by this format:
    options = {
        "":{
            "default":'',
            "type":'',
            "description": ""
        },
    }

    """
    @classmethod
    def is_valid_path(self, path):
        if os.path.exists(path):
            return True
        else:
            return False
        
    def get_all_options_types(self):
        """
        This needs to be a method that grabs attributes from each class, and figures out what type they are as a TypeSafe object
        Also grab their default values
        """

        return self.options

    def initialize_from_dict(self, dict: dict, name):
        for item, info in dict.items(): 
            setattr(self, item, TypeSafe(default=info["default"], 
                                         type=info["type"], 
                                         object_name=item))
        logger.info(f'All options for {name} have been initialized to default values.')

    def todict(self):
        dict = {}
        for item in self.options.keys():
            dict[item] = getattr(self, item)
        returndict = copy.deepcopy(dict)
        return returndict
    
    def __getattribute__(self, attribute):
        dataobject = object.__getattribute__(self, attribute)
        if isinstance(dataobject, TypeSafe):
            value = dataobject.value
        else:
            value = dataobject
        return value
    
    def __setattr__(self, attribute, value):
            if not hasattr(self, attribute):
                object.__setattr__(self, attribute, value)
            else:
                dataobject = object.__getattribute__(self, attribute)
                dataobject.__setattr__("value", value)
            self.checks()

class GeneralOptions(RTOptions):
    """Creates a class with general options relevant to any Ringtail process.
    """
    options = {
        "docking_mode":{
            "type": str,
            "default": "dlg",
            "description": 'specify AutoDock program used to generate results. Available options are "DLG" and "Vina". Vina mode will automatically change --pattern to *.pdbqt'
        },
        "print_summary":{
            "type": bool,
            "default": False,
            "description":"prints summary information about stored data to STDOUT."
        },  
        "logging_level":{
            "type": str,
            "default": "DEBUG",
            "description": '''
                            "WARNING": Prints errors and warnings to stout only. 
                            "INFO": Print results passing filtering criteria to STDOUT. NOTE: runtime may be slower option used
                            "DEBUG": Print additional error information to STDOUT
                            '''
        },  
    }
    
    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)
    
    def checks(self):
        if self.docking_mode not in ["dlg", "vina"]:
            logger.error(f'Docking mode {self.docking_mode} is not supported. Please choose between "dlg" and "vina".')

class InputFiles(RTOptions):
    """ Class that handles sources of data to be written including ligand data paths and how 
    to traverse them, and options to store receptor.
    """
    options = {
        "file":{
            "default":None,
            "type":list,
            "description": "ligand docking output file to save. Compressed (.gz) files allowed. Only 1 receptor allowed."
        },
        "file_path":{
            "default":None,
            "type":list,
            "description": " directory(s) containing docking output files to save. Compressed (.gz) files allowed"
        },
        "file_list":{
            "default":None,
            "type":list,
            "description": "file(s) containing the list of docking output files to save; relative or absolute paths are allowed. Compressed (.gz) files allowed"
        },
        "file_pattern":{
            "default":"*.dlg*",
            "type":str,
            "description": 'specify which pattern to use when searching for result files to process [only with "--file_path"]'
        },
        "recursive":{
            "default":None,
            "type":bool,
            "description": "enable recursive directory scan when --file_path is used"
        },
        "receptor_file":{
            "default":None,
            "type":str,
            "description": "Use with Vina mode. Give file for receptor PDBQT."
        },
        "save_receptor":{
            "default":None,
            "type":bool,
            "description": "Saves receptor PDBQT to database. Receptor location must be specied with in --file, --file_path directory or --file_list file"
        },
        "target":{
            "default":None,
            "type":str,
            "description": "name of receptor"
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        if hasattr(self, "target"): # ensures last item in the option dictionary has been
            if type(self.target) != str:
                if self.receptor_file is None:
                    pass
                elif self.receptor_file is not None and self.is_valid_path(self.receptor_file):
                    self.target = (os.path.basename(self.receptor_file).split(".")[0])
                else:
                    raise OptionError("The receptor PDBQT file path is not valid. Please check location of receptor file and --receptor_file option")
        
class ResultsProcessingOptions(RTOptions):
    """ Class that holds database write options that affects write time, such as how to 
    break up data files, number of computer processes to use, and and how many poses to store."""

    options = {
        "store_all_poses":{
            "default":False,
            "type":bool,
            "description": "Store all poses from input files. Overrides --max_poses"
        },
        "max_poses":{
            "default":3,
            "type":int,
            "description": "Store top pose for top n clusters"
        },
        "add_interactions":{
            "default":False,
            "type":bool,
            "description": "Find interactions between ligand poses and receptor and save to database. Requires receptor PDBQT to be given with input files (all modes) and --receptor_file to be specified with Vina mode. SIGNIFICANTLY INCREASES DATBASE WRITE TIME."
        },
        "interaction_tolerance":{
            "default":None,
            "type":float,
            "description": 'Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired. Only compatible with ADGPU mode'
        },
        "interaction_cutoffs":{
            "default":[3.7, 4.0],
            "type":list,
            "description": "Use with --add_interactions, specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. These are the default cutoffs."
        },
        "max_proc":{
            "default":None,
            "type":int,
            "description": "Maximum number of processes to create during parallel file parsing. Defaults to number of CPU processors."
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        # Check interaction_tolerance compatibilities 
        if hasattr(self, "max_proc"):
            if self.store_all_poses is True and self.interaction_tolerance is not None:
                logger.warning("Cannot use interaction_tolerance with store_all_poses. Removing interaction_tolerance.")
                self.interaction_tolerance = None       
        
class StorageOptions(RTOptions):
    """ Class that handles options for the storage (database) manager class, including
    conflict handling and result clustering and ordering.
    """

    options = {
        "append_results":{
            "default":None,
            "type":bool,
            "description": "Add new results to an existing database, specified by database choice in ringtail initialization or --input_db in cli."
        },
        "duplicate_handling":{
            "default":None,
            "type":str,
            "description": "specify how duplicate Results rows should be handled when inserting into database. Options are 'ignore' or 'replace'. Default behavior will allow duplicate entries."
        },
        "filter_bookmark":{
            "default":None,
            "type":str,
            "description": "Perform filtering over specified bookmark. (in output group in CLI)"
        },
        "overwrite":{
            "default":None,
            "type":bool,
            "description": "by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database"
        },
        "order_results":{
            "default":None,
            "type":str,
            "description": '''Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION.
                            available fields are:  
                            "e" (docking_score), 
                            "le" (ligand efficiency), 
                            "delta" (delta energy from best pose), 
                            "ref_rmsd" (RMSD to reference pose), 
                            "e_inter" (intermolecular energy), 
                            "e_vdw" (van der waals energy), 
                            "e_elec" (electrostatic energy), 
                            "e_intra" (intermolecular energy), 
                            "n_interact" (number of interactions), 
                            "rank" (rank of ligand pose), 
                            "run" (run number for ligand pose), 
                            "hb" (hydrogen bonds); '''
        },
        "outfields":{
            "default":"Ligand_name,e",
            "type":str,
            "description": '''defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are: \n
                            "Ligand_name" (Ligand name), 
                            "e" (docking_score), 
                            "le" (ligand efficiency), 
                            "delta" (delta energy from best pose), 
                            "ref_rmsd" (RMSD to reference pose), 
                            "e_inter" (intermolecular energy), 
                            "e_vdw" (van der waals energy), 
                            "e_elec" (electrostatic energy), 
                            "e_intra" (intermolecular energy), 
                            "n_interact" (number of iteractions), 
                            "ligand_smile" , 
                            "rank" (rank of ligand pose), 
                            "run" (run number for ligand pose), 
                            "hb" (hydrogen bonds), 
                            "receptor" (receptor name); '''
        },
        "output_all_poses":{
            "default":None,
            "type":bool,
            "description": "By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged."
        },
        "mfpt_cluster":{
            "default":None,
            "type":float,
            "description": "Cluster filtered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands."
        },
        "interaction_cluster":{
            "default":None,
            "type":float,
            "description": "Cluster filtered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions."
        },
        "bookmark_name":{
            "default":"passing_results",
            "type":str,
            "description": "name for resulting book mark file. Default value is 'passing_results'"
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        
        if hasattr(self, "bookmark_name"):
            # Makes sure all default values have been set once, so comparisons can start
            if self.duplicate_handling is not None:
                if self.duplicate_handling.upper() not in ["IGNORE", "REPLACE"]:
                    logger.warning(
                        f"--duplicate_handling option {self.duplicate_handling} not allowed. Reverting to default behavior."
                    )
                    self.duplicate_handling = None
            if self.order_results is not None and self.order_results not in self.order_options:
                raise OptionError(
                    "Requested ording option that is not available. Please see --help for available options."
                )
            # Make sure we include ligand name in output columnds
            if self.outfields is not None and "Ligand_name" not in self.outfields:
                self.outfields = "Ligand_name," + self.outfields
        
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

class ReadOptions(RTOptions):
    """ Class that holds options related to reading from the database, including format for
    result export and alternate ways of displaying the data (plotting),
    """
    options = {
        "filtering":{
            "default":None,
            "type":bool,
            "description": "switch for whether or not filtering is to be performed, to accommodate cmdline tools."
        },
        "plot":{
            "default":None,
            "type":bool,
            "description": "Makes scatterplot of LE vs Best Energy, saves as scatter.png."
        },
        "find_similar_ligands":{
            "default":None,
            "type":str,
            "description": "Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering."
        },
        "export_bookmark_csv":{
            "default":None,
            "type":str,
            "description": "Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables"
        },
        "export_bookmark_db":{
            "default":None,
            "type":bool,
            "description": "Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db"
        },
        "export_query_csv":{
            "default":None,
            "type":str,
            "description": "Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]"
        },
        "export_receptor":{
            "default":None,
            "type":bool,
            "description": "Export stored receptor pdbqt. Will write to current directory."
        },
        "data_from_bookmark":{
            "default":None,
            "type":bool,
            "description": "Write log of --outfields data for bookmark specified by --bookmark_name. Must use without any filters."
        },
        "pymol":{
            "default":None,
            "type":bool,
            "description": "Lauch PyMOL session and plot of ligand efficiency vs docking score for molecules in bookmark specified with --bookmark_name. Will display molecule in PyMOL when clicked on plot. Will also open receptor if given."
        },
        "enumerate_interaction_combs":{ #TODO does this belong in filters?
            "default":None,
            "type":bool,
            "description": "When used with `max_miss` > 0, will log ligands/poses passing each separate interaction filter combination as well as union of combinations. Can significantly increase runtime."
        },
        "log_file":{
            "default":"output_log.txt",
            "type":str,
            "description": "by default, results are saved in 'output_log.txt'; if this option is used, ligands and requested info passing the filters will be written to specified file"
        },
        "export_sdf_path":{
            "default":"",
            "type":str,
            "description": "specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified."
        },
    }

    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        if hasattr(self, "export_sdf_path"):
            if self.export_sdf_path is not None and not self.export_sdf_path == "" and not self.export_sdf_path.endswith("/"): 
                self.export_sdf_path += "/"

class Filters(RTOptions):
    """
    Object that holds all optional filters.
    """

    options = {
        "eworst":{
            "default":None,
            "type":float,
            "description": "specify the worst energy value accepted"
        },
        "ebest":{
            "default":None,
            "type":float,
            "description": "specify the best energy value accepted"
        },
        "leworst":{
            "default":None,
            "type":float,
            "description": "specify the worst ligand efficiency value accepted"
        },
        "lebest":{
            "default":None,
            "type":'',
            "description": "specify the best ligand efficiency value accepted"
        },
        "score_percentile":{
            "default":None,
            "type":float,
            "description": "specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent."
        },
        "le_percentile":{
            "default":None,
            "type":float,
            "description": "specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent."
        },
        "vdw_interactions":{
            "default":[],
            "type": list,
            "description": "define van der Waals interactions with residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]"
        },
        "hb_interactions":{
            "default":[],
            "type":list,
            "description": "define HB (ligand acceptor or donor) interaction as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]"
        },
        "reactive_interactions":{
            "default":[],
            "type":list,
            "description": "check if ligand reacted with specified residue as [-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]. E.g., [('A:VAL:279:', True), ('A:LYS:162:', True)] -> [('chain:resname:resid:atomname', <wanted (bool)>), ('chain:resname:resid:atomname', <wanted (bool)>)]"
        },
        "interactions_count":{
            "default":[],
            "type":list,
            "description": "accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions. E.g., [('hb_count', 5)]"
        },
        "react_any":{
            "default":None,
            "type":bool,
            "description": "check if ligand reacted with any residue"
        },
        "max_miss":{
            "default":0,
            "type":int,
            "description": "Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with 'enumerate_interaction_combs' for enumeration of poses passing each individual combination of interaction filters."
        },
        "ligand_name":{
            "default":[],
            "type":list,
            "description": "specify ligand name(s). Will combine name filters with OR"
        },
        "ligand_substruct":{
            "default":[],
            "type":list,
            "description": "SMARTS pattern(s) for substructure matching"
        },
        "ligand_substruct_pos":{
            "default":[],
            "type":list,
            "description": "SMARTS pattern(s) for substructure matching, e.g., [''[Oh]C' 0 1.2 -5.5 10.0 15.5'] -> ['smart_string index_of_positioned_atom cutoff_distance x y z']"
        },
        "ligand_max_atoms":{
            "default":None,
            "type":int,
            "description": "Maximum number of heavy atoms a ligand may have"
        },
        "ligand_operator":{
            "default":"OR",
            "type":str,
            "description": "logical join operator for multiple SMARTS (default: OR)"
        },
    }
    def __init__(self):
        super().initialize_from_dict(self.options, self.__class__.__name__)

    def checks(self):
        """
        Ensures all values are internally consistent and valid. Runs once after all values are set initially,
        then every time a value is changed.
        
        """
        if not hasattr(self, "ligand_operator"):
            # means the object is not fully initialized yet
            pass
        else:
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
    