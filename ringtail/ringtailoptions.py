#TODO ringtail should have validation checks of all parameters at the time they are set, 
# so they can set values independently 

#TODO the interface needs to host the logger

#TODO what values are needed for what specific job

#TODO so there are two things going on, one is the group during cmd inputs,
# the other is regrouping them for what managers the options go to

#TODO start setting some limits for various parameters, maybe an acronym attribute
# that can be used to join with cmd line options, as well as ensuring default values
# are set when they need to.
# Then I can make do with the default methods in RT core
#TODO create a new method in RTCore to display what the default values are, organized 
# by dictionary/grouping

class rtargs: 
    '''Base class for all options used in Ringtail. Each option is set as its
    own object, and has metadata (some which is legacy from cmd line parser).
    The value of each option is safely set according to assigned type.'''

    def __init__(self, type):
        self.name = None
        self.altname = None
        self.dest = None 
        self.nargs = None
        self.default = None
        self.type = type
        self.choices = None
        self.required = False
        self.help = None 
        self.metavar = None

    def _set_value(self, value):
        '''Checks that set value is of correct type'''
        if not isinstance(value, self.type):
            raise TypeError(self.name + " must be set to type " + str(self.type))
        self.__value = value
    
    def _get_value(self):
        '''Checks if value has been set'''
        try: 
            return self.__value
        except:
            raise NameError(self.name + " does not have a value.")
        
    value = property(_get_value, _set_value)

    def set_attr(self, **kwargs):
        '''assigns keyword arguments for each argument object, used e.g., for CLOPtionparser'''
        for k, v in kwargs.items():
            setattr(self, k, v)

class input_db(rtargs):
    def __init__(self):
        super(input_db, self).__init__(str)
        input_db.set_attr(self, **{"name": __class__.__name__, 
                                    "altname":"i", 
                                    "help":"specify a database file to perform actions with", 
                                    "action":"store",
                                    "metavar":"DATABASE"})
        
class output_db(rtargs):
    def __init__(self):
        super(input_db, self).__init__(str)
        input_db.set_attr(self, **{"name": __class__.__name__, 
                                   "altname": "o",
                                    "help":"Name for output database file",
                                    "action":"store",
                                    "type":str,
                                    "metavar":"[FILE_NAME].DB",
                                    "default":"output.db"})

"""
    mode  = "dlg" #rw, rman

    output_db = "output.db" #w
    input_db = None #rw
    input_db = types.SimpleNamespace
    kwargs = {"altname":"i", "help":"specify a database file to perform actions with", "action":"store","type":"str","metavar":"DATABASE"}
    

    bookmark_name = None #r

    summary = None #rw
    verbose = False #rw
    debug = False #rw

class RTWrite(RTArgs):  
    max_poses = None #rman also rw? 
    store_all_poses = None #rman
    interaction_tolerance = 0.8 #rman
    add_interactions = None #rman
    interaction_cutoffs = None #rman
    receptor_file = None #rman
    save_receptor = None #rman implicit
    max_proc = None #rman
    #I wonder of these file options could be internal properties
    file = None #rman
    file_path = [] #rman 
    file_list = None #rman
    file_pattern = "*.dlg*" #rman
    recursive = None #part of file sources
    # file handling options
    append_results = None #storage
    duplicate_handling = None #storage, part of conflict handling
    overwrite = None #storage

class RTRead(RTArgs):
    results_view_name = None
    # are read only according to CLOptionparser
    find_similar_ligands = None
    plot = None
    pymol = None
    export_bookmark_db = None
    export_receptor = None
    data_from_bookmark = None
    export_bookmark_csv = None
    export_query_csv = None
    log_file = None #outman
    enumerate_interaction_combs = None

class RTFilters(RTRead):
    # goes into all_filters
    van_der_waals = None #optional filter
    hydrogen_bond = None #optional filter
    reactive_res = None #optional filter
    hb_count = None
    react_any = None #boolean
    max_miss = None
    #ligand
    name = None #optional filter
    max_nr_atoms = None #optional filter
    smarts = None #optional filter
    smarts_idxyz = None #optional filter
    smarts_join = None
    #"properties group" --> goes into all_filters
    eworst = None #optional filter
    ebest = None #optional filter
    leworst = None #optional filter
    lebest = None #optional filter
    score_percentile = None #optional filter
    le_percentile = None #optional filter
    # output options
    outfields = None
    order_results = None
    output_all_poses = None
    mfpt_cluster = None
    interaction_cluster = None
    export_sdf_path = None #outman
    filter_bookmark = None


# storage_opts
'''         "db_file": dbFile,
            "order_results": parsed_opts.order_results,
            "outfields": parsed_opts.outfields,
            "filter_bookmark": parsed_opts.filter_bookmark,
            "output_all_poses": parsed_opts.output_all_poses,
            "mfpt_cluster": parsed_opts.mfpt_cluster,
            "interaction_cluster": parsed_opts.interaction_cluster,
            "results_view_name": parsed_opts.results_view_name,
            "overwrite": parsed_opts.overwrite,
            "append_results": parsed_opts.append_results,
            "conflict_opt": conflict_handling,
        }'''
    

#rman_opts
'''         "chunk_size": 1,
            "mode": parsed_opts.mode,
            "store_all_poses": parsed_opts.store_all_poses,
            "max_poses": parsed_opts.max_poses,
            "interaction_tolerance": parsed_opts.interaction_tolerance,
            "add_interactions": parsed_opts.add_interactions,
            "interaction_cutoffs": parsed_opts.interaction_cutoffs,
            "receptor_file": parsed_opts.receptor_file,
            "target": receptor,
            "file_sources": file_sources,
            "file_pattern": parsed_opts.file_pattern,
            "parser_manager": "multiprocessing",
            "max_proc": parsed_opts.max_proc,
        }'''

#outman_opts
'''{
            "log_file": parsed_opts.log_file,
            "export_sdf_path": parsed_opts.export_sdf_path,
        }'''


"""