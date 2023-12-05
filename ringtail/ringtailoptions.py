from ringtail import RingtailCore, APIOptionParser
#TODO ringtail should have validation checks of all parameters at the time they are set, 
# so they can set values independently 

class RingtailArguments(object): 
    # set all variables to default values
    process_mode = None
    # write args
    output_db = "output.db"
    input_db = None
    bookmark_name = None
    mode  = "dlg"
    summary = None
    verbose = None
    debug = None
    file = None
    file_path = []
    file_list = None
    file_pattern = "*.dlg*"
    recursive = None
    append_results = False
    duplicate_handling = None
    save_receptor = None
    overwrite = False
    max_poses = None
    store_all_poses = None
    interaction_tolerance = None
    add_interactions = None
    interaction_cutoffs = None
    receptor_file = None
    max_proc = None
    # Interaction args
    van_der_waals = None
    hydrogen_bond = None
    reactive_res = None
    hb_count = None
    react_any = None
    max_miss = None
    enumerate_interaction_combs = None
    # ligand_group
    name = None
    max_nr_atoms = None
    smarts = None
    smarts_idxyz = None
    smarts_join = None
    # properties_group
    eworst = None
    ebest = None
    leworst = None
    lebest = None
    score_percentile = None
    le_percentile = None
    # read_parser
    results_view_name = None
    summary = None
    verbose = None
    debug = None
    # output group
    log_file = None
    outfields = None
    order_results = None
    output_all_poses = None
    mfpt_cluster = None
    interaction_cluster = None
    export_bookmark_csv = None
    export_query_csv = None
    export_sdf_path = None
    export_bookmark_db = None
    export_receptor = None
    data_from_bookmark = None
    filter_bookmark = None
    find_similar_ligands = None
    plot = None
    pymol = None

    def addresults(self, rtcore):
        defaults = RingtailCore.get_defaults()
        with rtcore:
            rtcore.add_results()


    
# #TODO this needs to go in unit test instead
# rtopt = RingtailArguments()
# rtopt.process_mode = "write"
# rtopt.output_db = "newdb.db"
# rtopt.file_path = [['/Users/maylinnp/forlilab/Ringtail/test/test_data/']]
# rtopt.recursive = True

# rtcore = RingtailCore()
# APIopts = APIOptionParser(ringtail_core=rtcore, opts =rtopt)
# rtopt.addresults(rtcore)