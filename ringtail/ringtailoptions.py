#TODO the interface needs to host the logger

class RTArgs: 
    '''Base class for all options used in Ringtail. Each option is set as its
    own object, and has metadata (some which is legacy from cmd line parser).
    The value of each option is safely set according to assigned type.'''

    def __init__(self, type):
        self.type = type
        self.value = None

    @property 
    def value(self):
        '''Checks if value has been set'''
        try:
            return self._value
        except: 
            raise NameError(str(self.name) + " does not have a value.")  
   
    @value.setter
    def value(self, val):
        '''Checks that set value is of correct type'''
        if not isinstance(val, self.type) and val is not None:
            raise TypeError(str(self.name) + " must be set to type " + str(self.type))
        self._value = val

    def _set_attr(self, **kwargs):
        '''assigns keyword arguments for each argument object, used e.g., for CLOPtionparser'''
        for k, v in kwargs.items():
            setattr(self, k, v)

# # # Public classes and methods
def parserargs(object):
    """Parses the attributes of specified object to a string used to build the cmd line parser"""

    parser_string = '"-' + object.altname + '", "--' + object.name + '",'
    attrlist = [i for i in vars(object) if 'name' not in i] # remove attributes already assigned
    for var in attrlist:
        value = getattr(object, var)
        if value != None:
            parser_string += str(var) + '="' + str(value) + '",'
    return parser_string

class Filters():
    """ Class that holds all ringtail filters and their values"""
    def __init__(self, 
                 eworst = None,
                 ebest = None,
                 leworst = None,
                 lebest = None,
                 score_percentile = None,
                 le_percentile = None,
                 name = None,
                 max_nr_atoms = None,
                 smarts = None,
                 smarts_idxyz = None,
                 smarts_join = None,
                 van_der_waals = None,
                 hydrogen_bond = None,
                 reactive_res = None,
                 hb_count = None,
                 react_any = None,
                 max_miss = None,
                 enumerate_interaction_combs = None):
        self.eworst = EWorst(eworst)
        self.ebest = EBest(ebest)
        self.leworst = LEWorst(leworst),
        self.lebest = LEBest(lebest)
        self.score_percentile = ScorePercentile(score_percentile)
        self.le_percentile = LEPercentile(le_percentile)
        self.name = Name(name)
        self.max_nr_atoms = MaxNrAtoms(max_nr_atoms)
        self.smarts = Smarts(smarts)
        self.smarts_idxyz = SmartsIDXYZ(smarts_idxyz)
        self.smarts_join = SmartsJoin(smarts_join)
        self.van_der_waals = VanDerWaals(van_der_waals)
        self.hydrogen_bond = HydrogenBond(hydrogen_bond)
        self.reactive_res = ReactiveRes(reactive_res)
        self.hb_count = HBCount(hb_count)
        self.react_any = ReactAny(react_any)
        self.max_miss = MaxMiss(max_miss)
        self.enumerate_interaction_combs = EnumerateInteractionCombs(enumerate_interaction_combs)
    

# # # General options for Ringtail
class ProcessMode(RTArgs): #TODO this might need some finagling 
    def __init__(self, value=None):
        super(ProcessMode, self).__init__(str)
        ProcessMode._set_attr(self, **{"name":"process_mode",
                                    "help":'f"Specify if should write to or read from database. To show options of each mode use the in-line help, e.g.: {os.path.basename(__main__.__file__)} read -h"', 
                                    "dest":"process_mode"})
        if value != None:
            self.value = value

class Mode(RTArgs):
    def __init__(self, value=None):
        super(Mode, self).__init__(str)
        Mode._set_attr(self, **{"name":"mode", 
                                "altname":"m",
                                "help":'specify AutoDock program used to generate results. Available options are "DLG" and "Vina". Vina mode will automatically change --pattern to *.pdbqt',
                                "action":"store",
                                "type":str,
                                "metavar":"[dlg] or [vina]"})
        if value != None:
                self.value = value

class Summary(RTArgs):
    def __init__(self, value=None):
        super(Summary, self).__init__(bool)
        Summary._set_attr(self, **{"name":"summary", 
                                    "altname":"su",
                                    "help":'Prints summary information about stored data to STDOUT. Includes number of stored ligands and poses, min and max docking score and ligand efficiency, and 1%% (percentile) and 10%% (percentile) energy and ligand efficiency.',
                                    "action":"store_true"})
        if value != None:
            self.value = value
        
class Verbose(RTArgs):
    def __init__(self, value=None):
        super(Verbose, self).__init__(bool)
        Verbose._set_attr(self, **{"name":"verbose", 
                                    "altname":"v",
                                    "help":'Print results passing filtering criteria to STDOUT. NOTE: runtime may be slower option used.',
                                    "action":"store_true"})
        if value != None:
            self.value = value
        
class Debug(RTArgs):
    def __init__(self, value=None):
        super(Debug, self).__init__(bool)
        Debug._set_attr(self, **{"name":"debug", 
                                "altname":"d",
                                "help":'Print additional error information to STDOUT.',
                                "action":"store_true"})
        if value != None:
            self.value = value
       
# # # Database and programming classes
class InputDB(RTArgs): #can be source of data for write, and only source of data for read
    def __init__(self, value=None):
        super(InputDB, self).__init__(str)
        InputDB._set_attr(self, **{"name":"input_db", 
                                    "altname":"i", 
                                    "help":"specify a database file to perform actions with", 
                                    "action":"store",
                                    "metavar":"DATABASE"})
        if value != None:
            self.value = value
        
class OutputDB(RTArgs): #only used in write. I am wondering if there can be just one db name, and 
    # a switch on whether or not results should be appended
    def __init__(self, value=None):
        super(OutputDB, self).__init__(str)
        OutputDB._set_attr(self, **{"name":"output_db",
                                   "altname":"o",
                                    "help":"Name for output database file",
                                    "action":"store",
                                    "metavar":"[FILE_NAME].DB",
                                    "default":"output.db"})
        if value != None:
            self.value = value
    
class BookmarkName(RTArgs): 
    def __init__(self, value=None, process_mode="write"):
        super(BookmarkName, self).__init__(str)
        if process_mode == "read":
            BookmarkName._set_attr(self, **{"name":"bookmark_name", 
                                            "altname":"s",
                                            "help":'Specify name for db view of passing results to create or export from',
                                            "action":"store",
                                            "metavar":"STRING",
                                            "dest":"results_view_name"})
        else:
            BookmarkName._set_attr(self, **{"name":"bookmark_name",
                                            "altname":"s",
                                            "help":'Specify name for db view of passing results to create (write mode) or export from (read mode)',
                                            "action":"store",
                                            "metavar":"STRING"})
        if value != None:
            self.value = value
        
class AppendResults(RTArgs):
    def __init__(self, value=None):
        super(AppendResults, self).__init__(bool)
        AppendResults._set_attr(self, **{"name":"append_results", 
                                    "altname":"a",
                                    "help":"Add new results to an existing database, specified by --input_db",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class DuplicateHandling(RTArgs):
    def __init__(self, value=None):
        super(DuplicateHandling, self).__init__(str)
        DuplicateHandling._set_attr(self, **{"name":"duplicate_handling", 
                                    "altname":"dh",
                                    "help":'specify how duplicate Results rows should be handled when inserting into database. Options are "ignore" or "replace". Default behavior will allow duplicate entries.',
                                    "action":"store",
                                    "metavar":"'ignore' or 'replace'"})
        if value != None:
            self.value = value

class Overwrite(RTArgs):
    def __init__(self, value=None):
        super(Overwrite, self).__init__(bool)
        Overwrite._set_attr(self, **{"name":"overwrite", 
                                    "altname":"ov",
                                    "help":"by default, if a log file exists, it doesn't get overwritten and an error is returned; this option enable overwriting existing log files. Will also overwrite existing database",
                                    "action":"store_true"})
        if value != None:
            self.value = value
        
# # # File classes
class File(RTArgs):
    def __init__(self, value=None):
        super(File, self).__init__(list)
        File._set_attr(self, **{"name":"file", 
                                "altname":"f",
                                "help":"ligand docking output file to save. Compressed (.gz) files allowed. Only 1 receptor allowed.",
                                "action":"append",
                                "metavar":"FILENAME.[DLG/PDBQT][.gz]",
                                "nargs":"+"})
        if value != None:
            self.value = value
        
class FilePath(RTArgs):
    def __init__(self, value=None):
        super(FilePath, self).__init__(list)
        FilePath._set_attr(self, **{"name":"file_path", 
                                    "altname":"fp",
                                    "help":"directory(s) containing docking output files to save. Compressed (.gz) files allowed",
                                    "action":"append",
                                    "metavar":"DIRNAME",
                                    "nargs":"+"})
        if value != None:
            self.value = value

class FileList(RTArgs):
    def __init__(self, value=None):
        super(FileList, self).__init__(list)
        FileList._set_attr(self, **{"name":"file_list", 
                                    "altname":"fl",
                                    "help":"file(s) containing the list of docking output files to save; relative or absolute paths are allowed. Compressed (.gz) files allowed",
                                    "action":"append",
                                    "metavar":"FILENAME",
                                    "nargs":"+"})
        if value != None:
            self.value = value
        
class Pattern(RTArgs):
    def __init__(self, value=None):
        super(Pattern, self).__init__(str)
        Pattern._set_attr(self, **{"name":"pattern", 
                                    "altname":"p",
                                    "help":'specify which pattern to use when searching for result files to process [only with "--file_path"]',
                                    "action":"store",
                                    "metavar":"PATTERN",
                                    "dest":"file_pattern"})
        if value != None:
            self.value = value
        
class Recursive(RTArgs):
    def __init__(self, value=None):
        super(Recursive, self).__init__(bool)
        Recursive._set_attr(self, **{"name":"recursive", 
                                    "altname":"r",
                                    "help":"enable recursive directory scan when --file_path is used",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class SaveReceptor(RTArgs):
    def __init__(self, value=None):
        super(SaveReceptor, self).__init__(bool)
        SaveReceptor._set_attr(self, **{"name":"save_receptor", 
                                    "altname":"sr",
                                    "help":"Saves receptor PDBQT to database. Receptor location must be specied with '--receptor_file'",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class ReceptorFile(RTArgs):
    def __init__(self, value=None):
        super(ReceptorFile, self).__init__(str)
        ReceptorFile._set_attr(self, **{"name":"receptor_file", 
                                    "altname":"rf",
                                    "help":"Use with Vina mode. Give file for receptor PDBQT.",
                                    "action":"store",
                                    "metavar":"STRING"})
        if value != None:
            self.value = value

# # # Write options

class StoreAllPoses(RTArgs):
    def __init__(self, value=None):
        super(StoreAllPoses, self).__init__(bool)
        StoreAllPoses._set_attr(self, **{"name":"store_all_poses", 
                                    "altname":"sa",
                                    "help":"Store all poses from input files. Overrides --max_poses",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class MaxPoses(RTArgs):
    def __init__(self, value=None):
        super(MaxPoses, self).__init__(int)
        MaxPoses._set_attr(self, **{"name":"max_poses", 
                                    "altname":"mp",
                                    "help":"n: Store top pose for top n clusters",
                                    "action":"store",
                                    "metavar":"INT"})
        if value != None:
            self.value = value

class AddInteractions(RTArgs):
    def __init__(self, value=None):
        super(AddInteractions, self).__init__(bool)
        AddInteractions._set_attr(self, **{"name":"add_interactions", 
                                    "altname":"ai",
                                    "help":"Find interactions between ligand poses and receptor and save to database. Requires receptor PDBQT to be given with input files (all modes) and --receptor_file to be specified with Vina mode. SIGNIFICANTLY INCREASES DATBASE WRITE TIME.",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class InteractionTolerance(RTArgs):
    def __init__(self, value=None):
        super(InteractionTolerance, self).__init__(float)
        InteractionTolerance._set_attr(self, **{"name":"interaction_tolerance", 
                                    "altname":"it",
                                    "help":"Will add the interactions for poses within some tolerance RMSD range of the top pose in a cluster to that top pose. Can use as flag with default tolerance of 0.8, or give other value as desired. Only compatible with ADGPU mode",
                                    "action":"store",
                                    "metavar":"FLOAT",
                                    "const":0.8,
                                    "nargs":"?"})
        if value != None:
            self.value = value

class InteractionCutoffs(RTArgs):
    def __init__(self, value=None):
        super(InteractionCutoffs, self).__init__(list)
        InteractionCutoffs._set_attr(self, **{"name":"interaction_cutoffs", 
                                    "altname":"ic",
                                    "help":"Use with --add_interactions, specify distance cutoffs for measuring interactions between ligand and receptor in angstroms. Give as string, separating cutoffs for hydrogen bonds and VDW with comma (in that order). E.g. '-ic 3.7,4.0' will set the cutoff for hydrogen bonds to 3.7 angstroms and for VDW to 4.0. These are the default cutoffs.",
                                    "action":"store",
                                    "metavar":"[HB CUTOFF],[VDW CUTOFF]"})
        if value != None:
            self.value = value

class MaxProc(RTArgs):
    def __init__(self, value=None):
        super(MaxProc, self).__init__(int)
        MaxProc._set_attr(self, **{"name":"max_proc", 
                                    "altname":"mpr",
                                    "help":"Maximum number of processes to create during parallel file parsing. Defaults to number of CPU processors.",
                                    "action":"store"})
        if value != None:
            self.value = value

# Read options
# # # Output options, some filters but mostly file handling
#Output file setup
class LogFile(RTArgs):
    def __init__(self, value=None):
        super(LogFile, self).__init__(str)
        LogFile._set_attr(self, **{"name":"log_file", 
                                    "altname":"l",
                                    "help":'by default, results are saved in "output_log.txt"; if this option is used, ligands and requested info passing the filters will be written to specified file',
                                    "action":"store",
                                    "metavar":"[FILE_NAME].TXT"})
        if value != None:
            self.value = value

class Outfields(RTArgs):
    def __init__(self, value=None):
        super(Outfields, self).__init__(str)
        Outfields._set_attr(self, **{"name":"outfields", 
                                    "altname":"of",
                                    "help":
                                        'defines which fields are used when reporting the results (to stdout and to the log file); fields are specified as comma-separated values, e.g. "--outfields=e,le,hb"; by default, docking_score (energy) and ligand name are reported; ligand always reported in first column available fields are:  '
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
                                        "Fields are "
                                        "printed in the order in which they are provided. Ligand name will always be returned and will be added in first position if not specified.",
                                    "action":"store",
                                    "metavar":"FIELD1,FIELD2,..."})
        if value != None:
            self.value = value

class OrderResults(RTArgs):
    def __init__(self, value=None):
        super(OrderResults, self).__init__(str)
        OrderResults._set_attr(self, **{"name":"order_results", 
                                    "altname":"ord",
                                    "help":"Stipulates how to order the results when written to the log file. By default will be ordered by order results were added to the database. ONLY TAKES ONE OPTION."
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
                                        '"hb" (hydrogen bonds); ',
                                    "action":"store",
                                    "metavar":"STRING"})
        if value != None:
            self.value = value

#These are sort of filters
class OutputAllPoses(RTArgs):
    def __init__(self, value=None):
        super(OutputAllPoses, self).__init__(bool)
        OutputAllPoses._set_attr(self, **{"name":"output_all_poses", 
                                    "altname":"oap",
                                    "help":"By default, will output only top-scoring pose passing filters per ligand. This flag will cause each pose passing the filters to be logged.",
                                    "action":"store_true",})
        if value != None:
            self.value = value

class MFPTCLuster(RTArgs):
    def __init__(self, value=None):
        super(MFPTCLuster, self).__init__(float)
        MFPTCLuster._set_attr(self, **{"name":"mfpt_cluster", 
                                    "altname":"mfpc",
                                    "help":"Cluster filered ligands by Tanimoto distance of Morgan fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for selecting chemically dissimilar ligands.",
                                    "action":"store",
                                    "metavar":"FLOAT",
                                    "const": 0.5,
                                    "nargs":"?"})
        if value != None:
            self.value = value

class InteractionCluster(RTArgs):
    def __init__(self, value=None):
        super(InteractionCluster, self).__init__(float)
        InteractionCluster._set_attr(self, **{"name":"interaction_cluster", 
                                    "altname":"ifpc",
                                    "help":"Cluster filered ligands by Tanimoto distance of interaction fingerprints with Butina clustering and output ligand with lowest ligand efficiency from each cluster. Default clustering cutoff is 0.5. Useful for enhancing selection of ligands with diverse interactions.",
                                    "action":"store",
                                    "metavar":"FLOAT",
                                    "const": 0.5,
                                    "nargs":"?"})
        if value != None:
            self.value = value

class FindSimilarLigands(RTArgs):
    def __init__(self, value=None):
        super(FindSimilarLigands, self).__init__(str)
        FindSimilarLigands._set_attr(self, **{"name":"find_similar_ligands", 
                                    "altname":"fsl",
                                    "help":"Allows user to find similar ligands to given ligand name based on previously performed morgan fingerprint or interaction clustering.",
                                    "action":"store"})
        if value != None:
            self.value = value

#More output and export options
class ExportBookmarkCSV(RTArgs):
    def __init__(self, value=None):
        super(ExportBookmarkCSV, self).__init__(str)
        ExportBookmarkCSV._set_attr(self, **{"name":"export_bookmark_csv", 
                                    "altname":"xs",
                                    "help":"Create csv of the bookmark given with bookmark_name. Output as <bookmark_name>.csv. Can also export full database tables",
                                    "action":"store",
                                    "metavar":"BOOKMARK_NAME"})
        if value != None:
            self.value = value

class ExportQueryCSV(RTArgs):
    def __init__(self, value=None):
        super(ExportQueryCSV, self).__init__(str)
        ExportQueryCSV._set_attr(self, **{"name":"export_query_csv", 
                                    "altname":"xq",
                                    "help":"Create csv of the requested SQL query. Output as query.csv. MUST BE PRE-FORMATTED IN SQL SYNTAX e.g. SELECT [columns] FROM [table] WHERE [conditions]",
                                    "action":"store",
                                    "metavar":"[VALID SQL QUERY]"})
        if value != None:
            self.value = value

class ExportSDFPath(RTArgs):
    def __init__(self, value=None):
        super(ExportSDFPath, self).__init__(str)
        ExportSDFPath._set_attr(self, **{"name":"export_sdf_path", 
                                    "altname":"sdf",
                                    "help":"specify the path where to save poses of ligands passing the filters (SDF format); if the directory does not exist, it will be created; if it already exist, it will throw an error, unless the --overwrite is used  NOTE: the log file will be automatically saved in this path. Ligands will be stored as SDF files in the order specified.",
                                    "action":"store",
                                    "metavar":"DIRECTORY_NAME"})
        if value != None:
            self.value = value

class ExportBookmarkDB(RTArgs):
    def __init__(self, value=None):
        super(ExportBookmarkDB, self).__init__(bool)
        ExportBookmarkDB._set_attr(self, **{"name":"export_bookmark_db", 
                                    "altname":"xdb",
                                    "help":"Export a database containing only the results found in the bookmark specified by --bookmark_name. Will save as <input_db>_<bookmark_name>.db",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class ExportReceptor(RTArgs):
    def __init__(self, value=None):
        super(ExportReceptor, self).__init__(bool)
        ExportReceptor._set_attr(self, **{"name":"export_receptor", 
                                    "altname":"xr",
                                    "help":"Export stored receptor pdbqt. Will write to current directory.",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class DataFromBookmark(RTArgs):
    def __init__(self, value=None):
        super(DataFromBookmark, self).__init__(bool)
        DataFromBookmark._set_attr(self, **{"name":"data_from_bookmark", 
                                    "altname":"nd",
                                    "help":"Write log of --outfields data for bookmark specified by --bookmark_name. Must use without any filters.",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class FilterBookmark(RTArgs):
    def __init__(self, value=None):
        super(FilterBookmark, self).__init__(str)
        FilterBookmark._set_attr(self, **{"name":"filter_bookmark", 
                                    "altname":"fb",
                                    "help":"Perform filtering over specified bookmark.",
                                    "action":"store"})
        if value != None:
            self.value = value

class Plot(RTArgs):
    def __init__(self, value=None):
        super(Plot, self).__init__(bool)
        Plot._set_attr(self, **{"name":"plot", 
                                    "altname":"p",
                                    "help":"Makes scatterplot of LE vs Best Energy, saves as scatter.png.",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class PyMOL(RTArgs):
    def __init__(self, value=None):
        super(PyMOL, self).__init__(bool)
        PyMOL._set_attr(self, **{"name":"pymol", 
                                    "altname":"py",
                                    "help":"Lauch PyMOL session and plot of ligand efficiency vs docking score for molecules in bookmark specified with --bookmark_name. Will display molecule in PyMOL when clicked on plot. Will also open receptor if given.",
                                    "action":"store_true"})
        if value != None:
            self.value = value

# # # Filters
#Properties group, energies and efficiencies 
class EWorst(RTArgs):
    def __init__(self, value=None):
        super(EWorst, self).__init__(float)
        EWorst._set_attr(self, **{"name":"eworst", 
                                    "altname":"e",
                                    "help":"specify the worst energy value accepted",
                                    "action":"store",
                                    "metavar":"FLOAT"})
        if value != None:
            self.value = value

class EBest(RTArgs):
    def __init__(self, value=None):
        super(EBest, self).__init__(float)
        EBest._set_attr(self, **{"name":"ebest", 
                                    "altname":"eb",
                                    "help":"specify the best energy value accepted",
                                    "action":"store",
                                    "metavar":"FLOAT"})
        if value != None:
            self.value = value

class LEWorst(RTArgs):
    def __init__(self, value=None):
        super(LEWorst, self).__init__(float)
        LEWorst._set_attr(self, **{"name":"leworst", 
                                    "altname":"le",
                                    "help":"specify the worst ligand efficiency value accepted",
                                    "action":"store",
                                    "metavar":"FLOAT"})
        if value != None:
            self.value = value

class LEBest(RTArgs):
    def __init__(self, value=None):
        super(LEBest, self).__init__(float)
        LEBest._set_attr(self, **{"name":"lebest", 
                                    "altname":"leb",
                                    "help":"specify the best ligand efficiency value accepted",
                                    "action":"store",
                                    "metavar":"FLOAT"})
        if value != None:
            self.value = value

class ScorePercentile(RTArgs):
    def __init__(self, value=None):
        super(ScorePercentile, self).__init__(float)
        ScorePercentile._set_attr(self, **{"name":"score_percentile", 
                                    "altname":"pe",
                                    "help":"specify the worst energy percentile accepted. Express as percentage e.g. 1 for top 1 percent.",
                                    "action":"store",
                                    "metavar":"FLOAT"})
        if value != None:
            self.value = value

class LEPercentile(RTArgs):
    def __init__(self, value=None):
        super(LEPercentile, self).__init__(float)
        LEPercentile._set_attr(self, **{"name":"le_percentile", 
                                    "altname":"ple",
                                    "help":"specify the worst ligand efficiency percentile accepted. Express as percentage e.g. 1 for top 1 percent.",
                                    "action":"store",
                                    "metavar":"FLOAT"})
        if value != None:
            self.value = value

#Ligand filtesr
class Name(RTArgs):
    def __init__(self, value=None):
        super(Name, self).__init__(str)
        Name._set_attr(self, **{"name":"name", 
                                    "altname":"n",
                                    "help":"specify ligand name(s). Will combine name filters with OR",
                                    "action":"store",
                                    "metavar":"STRING",
                                    "nargs":"+"})
        if value != None:
            self.value = value

class MaxNrAtoms(RTArgs):
    def __init__(self, value=None):
        super(MaxNrAtoms, self).__init__(int)
        MaxNrAtoms._set_attr(self, **{"name":"max_nr_atoms", 
                                    "altname":"mna",
                                    "help":"Maximum number of heavy atoms a ligand may have",
                                    "action":"store",
                                    "metavar":"INT"})
        if value != None:
            self.value = value

class Smarts(RTArgs):
    def __init__(self, value=None):
        super(Smarts, self).__init__(str)
        Smarts._set_attr(self, **{"name":"smarts", 
                                    "altname":"smarts",
                                    "help":"SMARTS pattern(s) for substructure matching",
                                    "action":"store",
                                    "metavar":"STRING",
                                    "nargs":"+"})
        if value != None:
            self.value = value

class SmartsIDXYZ(RTArgs):
    def __init__(self, value=None):
        super(SmartsIDXYZ, self).__init__(str)
        SmartsIDXYZ._set_attr(self, **{"name":"smarts_idxyz", 
                                    "altname":"smarts_idxyz",
                                    "help":"SMARTS, index of atom in SMARTS, cutoff dist, and target XYZ coords",
                                    "action":"store",
                                    "metavar":"STRING",
                                    "nargs":"+"})
        if value != None:
            self.value = value

class SmartsJoin(RTArgs):
    def __init__(self, value=None):
        super(SmartsJoin, self).__init__(str)
        SmartsJoin._set_attr(self, **{"name":"smarts_join", 
                                    "altname":"sj",
                                    "choices":["AND", "OR"],
                                    "help":"logical operator for multiple SMARTS (default: OR)",
                                    "action":"store",
                                    "metavar":"STRING"})
        if value != None:
            self.value = value

#Interaction filters
class VanDerWaals(RTArgs):
    def __init__(self, value=None):
        super(VanDerWaals, self).__init__(str)
        VanDerWaals._set_attr(self, **{"name":"van_der_waals", 
                                    "altname":"vdw",
                                    "help":"define van der Waals interactions with residue",
                                    "action":"append",
                                    "metavar":"[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"})
        if value != None:
            self.value = value

class HydrogenBond(RTArgs):
    def __init__(self, value=None):
        super(HydrogenBond, self).__init__(str)
        HydrogenBond._set_attr(self, **{"name":"hydrogen_bond", 
                                    "altname":"hb",
                                    "help":"define HB (ligand acceptor or donor) interaction",
                                    "action":"append",
                                    "metavar":"[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"})
        if value != None:
            self.value = value

class ReactiveRes(RTArgs):
    def __init__(self, value=None):
        super(ReactiveRes, self).__init__(str)
        ReactiveRes._set_attr(self, **{"name":"reactive_res", 
                                    "altname":"r",
                                    "help":"check if ligand reacted with specified residue",
                                    "action":"append",
                                    "metavar":"[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"})
        if value != None:
            self.value = value

class HBCount(RTArgs):
    def __init__(self, value=None):
        super(HBCount, self).__init__(int)
        HBCount._set_attr(self, **{"name":"hb_count", 
                                    "altname":"hc",
                                    "help":"accept ligands with at least the requested number of HB interactions. If a negative number is provided, then accept ligands with no more than the requested number of interactions",
                                    "action":"store",
                                    "metavar":"NUMBER"})
        if value != None:
            self.value = value

class ReactAny(RTArgs):
    def __init__(self, value=None):
        super(ReactAny, self).__init__(bool)
        ReactAny._set_attr(self, **{"name":"react_any", 
                                    "altname":"ra",
                                    "help":"check if ligand reacted with any residue",
                                    "action":"store_true"})
        if value != None:
            self.value = value

class MaxMiss(RTArgs):
    def __init__(self, value=None):
        super(MaxMiss, self).__init__(int)
        MaxMiss._set_attr(self, **{"name":"max_miss", 
                                    "altname":"mm",
                                    "help":"Will compute all possible combinations of interaction filters excluding up to max_miss numer of interactions from given set. Default will only return union of poses interaction filter combinations. Use with --enumerate_interaction_combs for enumeration of poses passing each individual combination of interaction filters.",
                                    "action":"store",
                                    "metavar":"INTEGER"})
        if value != None:
            self.value = value

class EnumerateInteractionCombs(RTArgs):
    def __init__(self, value=None):
        super(EnumerateInteractionCombs, self).__init__(bool)
        EnumerateInteractionCombs._set_attr(self, **{"name":"enumerate_interaction_combs", 
                                    "altname":"eic",
                                    "help":"Use with max_miss. If used, will output ligands passing each individual combination of interaction filters with max_miss.",
                                    "action":"store_true"})
        if value != None:
            self.value = value         
