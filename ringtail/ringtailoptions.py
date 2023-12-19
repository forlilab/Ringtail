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

class DBOptions():

    def __init__(self, 
                input_db=None,
                output_db="output.db",
                bookmark_name=None,
                append_results=None, 
                duplicate_handling=None,
                overwrite=None):
        self.input_db = InputDB(input_db).value 
        self.output_db = OutputDB(output_db).value 
        self.bookmark_name = BookmarkName(bookmark_name, process_mode).value
        self.append_results = AppendResults(append_results).value
        self.duplicate_handling = DuplicateHandling(duplicate_handling).value
        self.overwrite = Overwrite(overwrite).value
    #outman_opts
'''{
            "log_file": parsed_opts.log_file,
            "export_sdf_path": parsed_opts.export_sdf_path,
        }'''


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
    def __init__(self, value=None, process_mode="read"):
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
        