import sys
import argparse
from glob import glob
import os
import fnmatch
import gzip
import numpy as np
import multiprocessing as mp
from bashplotlib.histogram import plot_hist
from matplotlib import pyplot as plt
import time
import parsers as ad_parser
import sqlite3
from collections import defaultdict
import subprocess
from functools import partial

######################################################################################

###################################
# help and Filter file parser #####
###################################
time0 = time.perf_counter()

description = """description something of something"""
usage = "READ THE MANUAL!"
name= "vs_results_display"
epilog="""

REQUIRED PACKAGES
        Requires numpy, multiprocessing, bashplotlib, matplotlib, sqlite3.

AUTHOR
        Written by Althea Hansel-Harris. Based on code by Stefano Forli, PhD and Andreas Tillack, PhD.

REPORTING BUGS
        Please report bugs to:
        AutoDock mailing list   http://autodock.scripps.edu/mailing_list

COPYRIGHT
        Copyright (C) 2021 Stefano Forli Laboratory, Center for Computational Structural Biology, 
                     The Scripps Research Institute.
        GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>
""" 

option_groups = [
        { "INPUT" : {
            'desc': 'Specify input data',
            'args': [ 
                ('--file', {
                    'help':'ligand DLG or DLG.gz (compressed) file to filter', 
                    'action':'append', 'type':str, 'metavar': "FILENAME.DLG[.gz]", 'required':False},
                    ),
                ('--file_path', {
                    'help': 'directory containing DLG or DLG.gz (compressed) files to filter', 
                    'action':'store', 'type':str, 'metavar':"DIRNAME"},
                    ),
                ('--file_list', {'help': 'file containing the list of DLG or DLG.gz (compressed) files to filter; relative or absolute paths are allowed', 
                    'action':'store', 'type':str, 'metavar':"FILENAME"},
                    ),
                ('--recursive', {'help': 'enable recursive directory scan when --dir is used', 
                    'action':'store_true', 'default': False },
                    ),
                ('--pattern', {'help': 'specify which pattern to use when searching for files to process in directories [only with "--dir", default "*.dlg*"]', 
                    'action':'store', 'type': str, 'metavar': 'PATTERN','default': "*.dlg*" },
                    ),
                ('--filters_file', {
                    'help':'specify a file containing the filter definitions; each line in the file must contain the definition of a single filter, using the keywords of each filter as variable names, e.g.: "eworst=-3.0", "vdw=A:THR:276". NOTE: properties defined here have the precedence on the other CLI options!', 
                    'action':'store', 'type':str, 'metavar': "FILTERS_FILE", 'default':None},
                 ),
                ('--sql_db', {
                    'help':'specify a file containing a sqlite database', 
                    'action':'store', 'type':str, 'metavar': "SQL_DB", 'default':None},
                 ),
                    ],
        },},

        {'OUTPUT': {
            'desc': 'Manage the type of data reported and where it is written.',
            'args': [ 
                ('--output_sql', {
                    'help':('Name for output SQLite file'), 
                    'action':'store', 'type':str, 'metavar': "OUTLIST.TXT", 'default':"output.db"},
                    ),
                ('--num_clusters', {
                    'help':('n: Store top pose for top n clusters'), 
                    'action':'store', 'type':int, 'default':3},
                    ),
                ('--log', {
                    'help':('by default, results are printed in the terminal (STDOUT); '
                        'if this option is used, ligands passing the filters will be written '
                        'to this file'), 
                    'action':'store', 'type':str, 'metavar': "OUTLIST.TXT", 'default':None},
                    ),
                ('--out_fields', {
                    'help':'defines which fields are used when reporting '
                    'the results (either to stdout or to the log file); in the output, values '
                    'are separated by a "|" symbol; fields are specified ' 
                    'as comma-separated values, e.g. "--out_fields=e,le,fname"; by '
                    'default, energies_binding (energy) and file name (fname) are reported; fname always reported in first column'
                    'available fields are:  '
                    '"e" (energies_binding), '
                    '"le" (ligand efficiency), '
                    '"delta" (delta energy from best pose), '
                    '"ref_rmsd" (RMSD to reference pose), '
                    '"e_inter" (intermolecular energy), '
                    '"e_vdw" (van der waals energy), '
                    '"e_elec" (electrostatic energy), '
                    '"e_intra" (intermolecular energy), '
                    '"n_interact" (number of interactions), '
                    '"interactions" (all interactions), '
                    '"fname" (full path filename), '
                    '"ligand_smile" , '
                    '"rank" (rank of ligand pose), '
                    '"run" (run number for ligand pose), '
                    #'"mname" (molecule name from the \'move\' kw in the DLG); '
                    #'"name" (basename of the dlg[.gz] file); '
                    '"hb" (hydrogen bonds); '
                    #'"hba" (hydrogen bond, ligand acceptor); '
                    #'"hbd" (hydrogen bond, ligand donor); '
                    #'"r" (reacted: 1 (yes), 0 (no));'
                    #'"ri" (reacted residue info); '
                    '"all" (print all poses passing filters). Fields are '
                    'printed in the order in which they are provided.',
                    'action':'store', 'type':str, 'metavar': "FIELD1,FIELD2,...", 
                    'default':'fname,e'},
                    ),
                ('--export_poses_path', {
                    'help':('specify the path where to save poses of ligands passing the filters (PDBQT format); '
                        'if the directory does not exist, it will be created; if it already exist, it will throw '
                        'an error, unless the --overwrite is used  [QUESTION????: NOTE: the log file will be automatically saved in this path.]'),
                    'action':'store', 'type':str, 'metavar': "DIRECTORY_NAME", 'default':None},
                    ),
                ('--no_header', {
                    'help':('by default, a commented header ("#") with a summary of the filters used '
                        'is written to both STDOUT and the log file; this option suppresses the header'),
                    'action':'store_true', 'default':False },
                    ),
                ('--no_print', {
                    'help':('suppress printing the results to STDOUT. NOTE: runtime may be faster if no_print option used.'),
                    'action':'store_true', 'default':False },
                    ),
                ('--plot', {
                    'help':('Makes best energy and le histograms and displays in command line. Makes scatterplot of LE vs Best Energy, saves as [filters_file].png or out.png if no filters_file given.'),
                    'action':'store_true', 'default':False},
                    ),
                ('--overwrite', {
                    'help':('by default, if a log file exists, it doesn\'t get '
                        'overwritten and an error is returned; this option enable overwriting existing log files'),
                    'action':'store_true', 'default':False },
                    ),


                ],
            },},

        {'PROPERTY FILTERS': {
            'desc': ('Specify energy and ligand efficiency filters'),
            'args': [ 
                ('--eworst', {
                    'help':'specify the worst energy value accepted', 
                    'action':'store', 'type':float, 'metavar': "FLOAT", 'default':-3.0},
                    ),
                ('--ebest', {
                    'help':'specify the best energy value accepted', 
                    'action':'store', 'type':float, 'metavar': "FLOAT",},
                    ),
                ('--leworst', {
                    'help':'specify the worst ligand efficiency value accepted', 
                    'action':'store', 'type':float, 'metavar': "FLOAT"},
                    ),
                ('--lebest', {
                    'help':'specify the best ligand efficiency value accepted', 
                    'action':'store', 'type':float, 'metavar': "FLOAT",},
                    )

                    ],}},

        {'LIGAND FILTERS': {
            'desc': ('Specify filters on ligands, including substructures or names'),
            'args': [ 
                ('--name', {
                    'help':'specify ligand name(s). Will seach OR', 
                    'action':'store', 'type':str, 'default':None},
                    ),
                ('--substructure', {
                    'help':'specify SMILES substring(s) to search for', 
                    'action':'store', 'type':str, 'default':None},
                    ),
                ('--substruct_flag', {
                    'help':'specify whether to search AND or OR for substructures. Default OR', 
                    'action':'store', 'type':str, 'default':"OR"},
                    )
                    ],}},
        {'INTERACTION FILTERS': {
            'desc': ('Specify interaction filters, either by count or by specific residue interaction. '
                     'Residue specifications are described using CHAIN:RES:NUM:ATOM_NAME, '
                     'and any combination is allowed, e.g.: CHAIN:::, :RES::, ::NUM:, :::ATOM_NAME, :RES:NUM:, etc... '
                     'Unwanted interactions can be defined by prepending "-" to the residue specification, e.g. "-B:THR:276:". '
                     'Multiple residues can be specified in a single option by separating them with a comma (e.g.: --vdw=B:THR:276:,B:HIS:226:), '
                     'or by repeating the interaction options (e.g. --vdw=B:THR:276: --vdw=B:HIS:226: ).'),
            'args': [ 
                ('--vdw', {
                    'help':'define van der Waals interactions with residue', 
                    'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                    ),
                ('--hb', {
                    'help':'define HB (ligand acceptor or donor) interaction',
                    'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                    ),
                ('--hba', {
                    'help':'define HB acceptor (ligand as acceptor) interaction',
                    'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                    ),
                ('--hbd', {
                    'help':'define HB donor (ligand as donor) interaction',
                    'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                    ),
                ('--hb_count', {
                    'help':('accept ligands with at least the requested number of HB interactions. '
                        'If a negative number is provided, then accept ligands with no more than the '
                        ' requested number of interactions'),
                    'action':'store', 'type':int, 'metavar': "NUMBER"},
                    ),
                # ('--react_any', {
                #     'help':'check if ligand reacted with any residue',
                #     'action':'store_true'},
                #     ),
                ('--react_res', {
                    'help':'check if ligand reacted with specified residue',
                    'action':'append', 'type':str, 'metavar': "[-][CHAIN]:[RES]:[NUM]:[ATOM_NAME]"},
                    ),
                ('--react_count', {
                    'help':('accept ligands that reacted with at least the requested number of residues; '
                        'if a negative number is provided, then accept ligands that do not react '
                        'with more than the specified number of residues'),
                    'action':'store', 'type':int, 'metavar': "NUMBER"},
                    ),
                ],
            },},
    ]

def read_filter_file(fname):
        """ parse the filter file to define filters """
        opts = []
        with open(fname, 'r') as fp:
            for l in fp.readlines():
                l = l.strip()
                # remove empty lines
                if not len(l):
                    continue
                # remove comments
                if l[0] == "#":
                    continue
                opts.append("--%s" % l.strip())
        return opts

def process_options(parsed_opts):
    """ convert command line options to the dict of filters """
    # check that required input options are provided
    file_sources = {} #'file':None, 'files_path':None, 'file_list':None}
    file_sources['file'] = parsed_opts.file
    if not parsed_opts.file_path is None:
        file_sources['file_path'] = {'path':parsed_opts.file_path, 'pattern':parsed_opts.pattern, 'recursive': parsed_opts.recursive}
    else:
        file_sources['file_path'] = None
    file_sources['file_list'] = parsed_opts.file_list
    print("")
    if (file_sources['file'] is None) and (file_sources['file_path'] is None) and (file_sources['file_list'] is None) and (parsed_opts.sql_db is None):
        print("*ERROR* at least one input option needs to be used:  --file, --file_path, --file_list, --sql_db")
        sys.exit(1)
    ##### parse output options
    output = {'fields': parsed_opts.out_fields,
              'log': parsed_opts.log,
              'header': not parsed_opts.no_header,
              'stdout': not parsed_opts.no_print,
              'overwrite': parsed_opts.overwrite,
              'export_poses_path': parsed_opts.export_poses_path,
            }
    # if a path for saving poses is specified, then the log will be written there
    if not parsed_opts.export_poses_path is None:
        output['log'] = "%s%s%s" % (output['export_poses_path'], os.sep, output['log'])
    #print(">>>>>>>>>>>>>chekc that fields are recognized")
    if (parsed_opts.log is None) and (parsed_opts.no_print is True):
        print("*ERROR* print to STDOUT is disabled and no log file has been specified; at least one output source needs to be used.")
        sys.exit(1)
    ##### filters
    # property filters
    properties = {'eworst':None, 'ebest':None, 'leworst':None, 'lebest':None}
    for kw, _ in properties.items():
        properties[kw] = getattr(parsed_opts, kw)
    # interaction filters (residues)
    interactions = {}
    res_interactions_kw = [('vdw', 'V'), ('hb','H'), ('hba','HA'), ('hbd','HD'), ('react_res','R')]
    for opt, _type in res_interactions_kw:
        interactions[_type] = []
        res_list = getattr(parsed_opts, opt)
        if res_list is None:
            continue
        found_res = []
        for res in res_list:
            if "," in res:
                for r in res.split(','):
                    found_res.append(r)
            else:
                found_res.append(res)
        for res in found_res:
            wanted = True
            if not res.count(":") == 3:
                print( ('*ERROR* [%s]: to specify a residue use '
                    'the format CHAIN:RES:NUM:ATOM_NAME. Any item can be omitted, '
                    'as far as the number of semicolons is always 3 '
                    '(e.g.: CHAIN:::, :RES::, CHAIN::NUM:, etc.)') % res )
                sys.exit()
            if res[0] == '-':
                res=res[1:]
                wanted = False
            interactions[_type].append( (res, wanted) )
    # count interactions
    interactions_count = []
    count_kw = [('hb_count', ("hb_count")), ('react_count',('R'))]
    for kw, pool in count_kw:
        c = getattr(parsed_opts, kw, None)
        if c is None:
            continue
        interactions_count.append((pool, c))
    #make dictionary for ligand filters
    ligand_filters_kw = [('name', 'N'), ('substructure', 'S'), ('substruct_flag', 'F')]
    ligand_filters = {}
    filter_ligands_flag = True
    for kw, _type in ligand_filters_kw:
        ligand_filters[_type] = []
        ligand_filter_list = getattr(parsed_opts, kw)
        if ligand_filter_list is None:
            filter_ligands_flag = False
        else:
            ligand_filter_list = ligand_filter_list.split(",")
            for fil in ligand_filter_list:
                ligand_filters[_type].append(fil)

    # TODO combine everything into a single parameters dictionary to be passed to DockingResultManager (with all the kws in its init)
    filters = {'properties':properties, 'interactions':interactions, 'interactions_count':interactions_count, "ligand_filters":ligand_filters, 'filter_ligands_flag':filter_ligands_flag, 'max_miss':0} 
    return file_sources, filters, output

##################################################
############ database helper functions ###########
##################################################

def create_connection(db_file):
    con = None
    try:
        con = sqlite3.connect(db_file)
    except Exception as e:
        print(e)
    return con

def create_results_table(conn):
    sql_results_table = """CREATE TABLE Results (
        Pose_ID             INTEGER PRIMARY KEY AUTOINCREMENT,
        LigName             VARCHAR NOT NULL,
        ligand_smile        VARCHAR[],
        pose_rank           INT[],
        run_number          INT[],
        energies_binding    FLOAT(4),          -- 8 x 4 Bytes = 32 Bytes
        leff                FLOAT(4),
        deltas              FLOAT(4),
        cluster_rmsd        FLOAT(4),
        reference_rmsd      FLOAT(4),
        energies_inter      FLOAT(4),
        energies_vdw        FLOAT(4),
        energies_electro    FLOAT(4),
        energies_flexLig    FLOAT(4),
        energies_flexLR     FLOAT(4),
        energies_intra      FLOAT(4),
        energies_torsional  FLOAT(4),
        unbound_energy      FLOAT(4),
        nr_interactions     INT[],          -- 2 Bytes                                      (-> 46 Bytes)
        num_hb              INT[],
        about_x             FLOAT(4),          -- 3 x 4 Bytes = 12 Bytes
        about_y             FLOAT(4),
        about_z             FLOAT(4),
        trans_x             FLOAT(4),          -- 3 x 4 Bytes = 12 Bytes
        trans_y             FLOAT(4),
        trans_z             FLOAT(4),
        axisangle_x         FLOAT(4),          -- 4 x 4 Bytes = 16 Bytes                       (-> 30 Bytes)
        axisangle_y         FLOAT(4),
        axisangle_z         FLOAT(4),
        axisangle_w         FLOAT(4),
        dihedrals           VARCHAR[]          -- 8-ish Bytes per dihedral value
    );
    -- total number of Bytes estimated per row for average 50 interactions and 10 torsions:
    --     (64 + 76 + 50 * 21 + 10 * 8) Bytes = 1,270 Bytes
    -- for 50 ga_runs per result: (1,270 * 50) Bytes = 63,500 Bytes (62 KB)
    -- for 10,000 results (typical package size) this amounts to about 62 KB * 10,000 = 591 MB"""

    create_index = "CREATE INDEX index_ligand_run on Results (LigName, run_number)"


    try:
        cur = conn.cursor()
        cur.execute(sql_results_table)
        cur.execute(create_index)
    except Exception as e:
        print(e)

    cur.close()
    conn.close()

def create_ligands_table(conn):
    ligand_table = """CREATE TABLE Ligands (
        LigName             VARCHAR NOT NULL PRIMARY KEY,
        ligand_smile        VARCHAR NOT NULL,
        input_pdbqt         VARCHAR[],
        best_binding        FLOAT(4),
        best_run            SMALLINT[])"""

    try:
        cur = conn.cursor()
        cur.execute(ligand_table)
    except Exception as e:
        print(e)

    cur.close()
    conn.close()

def create_interaction_index_table(conn):
    """create table of data about each unique interaction"""
    interaction_index_table = """CREATE TABLE Interaction_indices (
        interaction_id      INTEGER PRIMARY KEY AUTOINCREMENT,
        interaction_type    VARCHAR[],
        rec_chain           VARCHAR[],
        rec_resname         VARCHAR[],
        rec_resid           VARCHAR[],
        rec_atom            VARCHAR[],
        rec_atomid          VARCHAR[])"""

    try:
        cur = conn.cursor()
        cur.execute(interaction_index_table)
    except Exception as e:
        print(e)

    cur.close()
    conn.close()

def create_interaction_bv_table(conn, num_discrete_interactions):
    interaction_columns = range(num_discrete_interactions)
    interact_columns_str = ""
    for i in interaction_columns:
        interact_columns_str += "Interaction_"+ str(i+1) + " INTEGER,\n"

    interact_columns_str = interact_columns_str.rstrip(",\n")

    bv_table = """CREATE TABLE Interaction_bitvectors (
    Pose_ID INTEGER PRIMARY KEY,
    {columns})""".format(columns = interact_columns_str)

    try:
        cur = conn.cursor()
        cur.execute(bv_table)
    except Exception as e:
        print(e)

    cur.close()
    conn.close()

def insert_results(conn, results_array):
    """takes array of database rows to insert, adds data to results table
    Inputs:
    -conn: database connection
    -ligand_array: numpy array of arrays"""

    print("Inserting results...")
    sql_insert = """INSERT INTO Results (
    LigName,
    ligand_smile,
    pose_rank,
    run_number,
    cluster_rmsd,
    reference_rmsd,
    energies_binding,
    leff,
    deltas,
    energies_inter,
    energies_vdw,
    energies_electro,
    energies_flexLig,
    energies_flexLR,
    energies_intra,
    energies_torsional,
    unbound_energy,
    nr_interactions,
    num_hb,
    about_x,
    about_y,
    about_z,
    trans_x,
    trans_y,
    trans_z,
    axisangle_x,
    axisangle_y,
    axisangle_z,
    axisangle_w,
    dihedrals
    ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"""

    try:
        cur = conn.cursor()
        cur.executemany(sql_insert, results_array)
        conn.commit()
        cur.close()
    except Exception as e:
        print("ERROR:")
        print(e)

    conn.close()

def insert_ligands(conn, ligand_array):
    print("Inserting ligand data...")
    sql_insert = '''INSERT INTO Ligands (
    LigName,
    ligand_smile,
    input_pdbqt,
    best_binding,
    best_run
    ) VALUES
    (?,?,?,?,?)'''

    try:
        cur = conn.cursor()
        cur.executemany(sql_insert, ligand_array)
        conn.commit()
        cur.close()
    except Exception as e:
        print(e)

    conn.close()

def insert_all_interactions(conn, all_interactions):
    print("Inserting interactions...")
    sql_insert = '''INSERT INTO Interaction_indices (
    interaction_type,
    rec_chain,
    rec_resname,
    rec_resid,
    rec_atom,
    rec_atomid
    ) VALUES (?,?,?,?,?,?)'''

    try:
        cur = conn.cursor()
        cur.executemany(sql_insert, all_interactions)
        conn.commit()
        cur.close()
    except Exception as e:
        print(e)

    conn.close()

def insert_interaction_BVs(conn, bitvectors, num_discrete_interactions):
    print("Inserting interaction bitvectors...")
    interaction_columns = range(num_discrete_interactions)
    column_str = ""
    filler_str = "?,"
    for i in interaction_columns:
        column_str += "Interaction_"+ str(i+1) + ", "
        filler_str += "?,"
    column_str = column_str.rstrip(", ")
    filler_str = filler_str.rstrip(",")

    sql_insert = '''INSERT INTO Interaction_bitvectors (Pose_ID, {columns}) VALUES ({fillers})'''.format(columns = column_str, fillers = filler_str)

    try:
        cur = conn.cursor()
        cur.executemany(sql_insert, bitvectors)
        conn.commit()
        cur.close()
    except Exception as e:
        print(e)

    conn.close()

def write_results_interaction_rows(ligand_dict):
    """writes list of lists of ligand values to be inserted into sqlite database"""

    ligand_rows = []
    interaction_rows = []

    #initialize list for sql row with name and smile string
    ligand_name = ligand_dict["ligname"]
    ligand_smile = ligand_dict["ligand_smile_string"]

    #print("Writing database rows for", ligand_name)

    ligand_data_keys = ["cluster_rmsds",
    "ref_rmsds",
    "scores",
    "leff",
    "delta",
    "intermolecular_energy",
    "vdw_hb_desolv",
    "electrostatics",
    "flex_ligand",
    "flexLigand_flexReceptor",
    "internal_energy",
    "torsional_energy",
    "unbound_energy"]

    ligand_interaction_keys = ["type",
    "chain",
    "residue",
    "resid",
    "recname",
    "recid"]

    stateVar_keys = ["pose_about",
    "pose_translations",
    "pose_quarternions"]
    
    ######get pose-specific data
    for i in range(len(ligand_dict)):
        #check if run is best for a cluster. We are only saving the top pose for each cluster
        pose_rank = i
        run_number = ligand_dict["sorted_runs"][pose_rank]
        try:
            cluster_top_pose_runs = ligand_dict["cluster_top_poses"][:parsed_opts.num_clusters] #will only select top n clusters. Default 3
        except IndexError:
            cluster_top_pose_runs = ligand_dict["cluster_top_poses"] #catch indexerror if not enough clusters for given ligand
        if run_number in cluster_top_pose_runs:
            ligand_data_list = [ligand_name, ligand_smile, pose_rank+1, run_number]
            #get energy data
            for key in ligand_data_keys:
                ligand_data_list.append(ligand_dict[key][pose_rank])

            #get interaction data
            """for key in interaction_keys:
                item = ligand_dict[key][pose_rank]
                if type(item) == list:
                    item_string = ""
                    for element in item:
                        item_string += element
                        item = item_string
                ligand_data_list.append(item)"""

            pose_interactions_dict = ligand_dict["interactions"][pose_rank]
            num_interactions = pose_interactions_dict["count"][0]
            ligand_data_list.append(num_interactions) 
            interaction_strings_list = []
            hb_flag = False
            for key in ligand_interaction_keys:
                interaction_data = pose_interactions_dict[key]
                if type(interaction_data) == list:
                    interaction_data_string = ""
                    for interaction in interaction_data:
                        interaction_data_string = interaction_data_string + interaction + ", "
                interaction_strings_list.append(interaction_data_string)
                if key == 'type' and not hb_flag: #add hb_count only once
                    hb_count = interaction_data_string.count("H")
                    ligand_data_list.append(hb_count)
                    hb_flag = True
            hb_flag = False
            #reformat interaction data to have all information for each interaction together
            interactions_string = ""
            for i in range(int(num_interactions)):
                single_interaction_string = ""
                for line in interaction_strings_list:
                    line_list = line.split(",")
                    single_interaction_string += line_list[i] + ":"
                interactions_string += single_interaction_string.replace(" ", "") + ", "
            interaction_rows.append(interactions_string)


            for key in stateVar_keys:
                stateVar_data = ligand_dict[key][pose_rank]
                for dim in stateVar_data:
                    ligand_data_list.append(dim)
            pose_dihedrals = ligand_dict["pose_dihedrals"][pose_rank]
            dihedral_string = ""
            for dihedral in pose_dihedrals:
                dihedral_string = dihedral_string + str(dihedral) + ", "
            ligand_data_list.append(dihedral_string)

            ligand_rows.append(ligand_data_list)

    return ligand_rows, interaction_rows

def write_ligand_row(ligand_dict):
    """writes row to be inserted into ligand table"""
    ligand_name = ligand_dict["ligname"]
    ligand_smile = ligand_dict["ligand_smile_string"]
    input_pdbqt = "\n".join(ligand_dict["ligand_input_pdbqt"])
    best_binding = ligand_dict["scores"][0]
    best_run = ligand_dict["sorted_runs"][0]

    return [ligand_name, ligand_smile, input_pdbqt, best_binding, best_run]

def get_all_interactions(interactions_strings, all_interactions_list, all_interactions_split):
    all_interactions_dict = {}
    for pose in interactions_strings:
        pose = pose.split(", ")
        for interaction in pose:
            if interaction.endswith(":"):
                    interaction = interaction.rstrip(":") #remove trailing colon
            if interaction not in all_interactions_dict and interaction != "":
                interaction_attributes= interaction.split(":")
                all_interactions_dict[interaction] = interaction
                all_interactions_list.append(interaction)
                all_interactions_split.append(interaction_attributes) #type,ligand atom, ligand id, chain, residue, resid, residue atom, residue atom id

    return all_interactions_list, all_interactions_split

def write_interaction_bitvector(pose_id, interactions_string, all_interactions):
    """takes string of interactions and all possible interactions and makes bitvector"""
    pose_bitvector = [pose_id+1]
    for interaction in all_interactions:
        if interaction in interactions_string:
            pose_bitvector.append(1) #true
        else:
            pose_bitvector.append(0) #false

    return pose_bitvector

def write_result_filtering_sql(filter_list):
    """ takes list of filters, writes sql filtering string"""

    #parse requested output fields and convert to column names in database
    output_fields = parsed_opts.out_fields.split(",")

    outfield_dict = {"e":"energies_binding",
                    "le":"leff",
                    "delta":"deltas",
                    "ref_rmsd":"reference_rmsd",
                    "e_inter":"energies_inter",
                    "e_vdw":"energies_vdw",
                    "e_elec":"energies_electro",
                    "e_intra":"energies_intra",
                    "n_interact":"nr_interactions",
                    "interactions":"interactions",
                    "fname":"LigName",
                    "ligand_smile":"ligand_smile",
                    "rank":"pose_rank",
                    "run":"run_number",
                    "hb":"num_hb"}
    outfield_string = ""
    for field in output_fields:
        if field in outfield_dict.keys():
            outfield_string += outfield_dict[field] + ", "
    outfield_string = outfield_string.rstrip(", ")

    interaction_filters = []

    sql_string = """SELECT {out_columns} FROM Results WHERE """.format(out_columns = outfield_string)

    for filter_tuple in filter_list:
        filter_key = filter_tuple[0]
        filter_value = filter_tuple[1]

        #write energy filters
        if filter_key == 'eworst':
            eworst_str = "energies_binding < {value} AND ".format(value = filter_value)
            sql_string += eworst_str

        if filter_key == 'ebest':
            ebest_str = "energies_binding > {value} AND ".format(value = filter_value)
            sql_string += ebest_str

        if filter_key == 'leworst':
            leworst_str = "leff < {value} AND ".format(value = filter_value)
            sql_string += leworst_str

        if filter_key == 'lebest':
            lebest_str = "leff > {value} AND ".format(value = filter_value)
            sql_string += lebest_str

        #write interaction filters
        if filter_key == 'hb_count':
            if filter_value > 0:
                hb_count_str = "hb_count > {value} AND ".format(value = filter_value)
                sql_string += hb_count_str
            else:
                hb_count_str = "hb_count < {value} AND ".format(value = -1*filter_value)
                sql_string += hb_count_str

        interaction_filter_keys = ["V",
        "H",
        "R"]
        #compile list of interactions to search for
        for key in interaction_filter_keys:
            if filter_key == key:
                for interact in filter_value:
                    interaction_string = key + ":" + interact[0]
                    interaction_filters.append(interaction_string.split(":"))

    interaction_filter_indices = []
    #for each interaction, get the index from the interactions_indices table
    for interaction in interaction_filters:
        interact_index_str = write_interaction_index_filtering_str(interaction)
        conn = create_connection(parsed_opts.output_sql)
        interaction_indices = select_from_db(conn, interact_index_str)
        for i in interaction_indices:
            interaction_filter_indices.append(i[0])

    #find pose ids for ligands with desired interactions
    if interaction_filter_indices != []:
        interaction_filter_str = write_interaction_filtering_str(interaction_filter_indices)
        conn = create_connection(parsed_opts.output_sql)
        interaction_pass_pose_ids = select_from_db(conn, interaction_filter_str)

        sql_string += "("
        for passed_pose in interaction_pass_pose_ids:
            passed_pose_str = "Pose_ID = {pose} OR ".format(pose = passed_pose[0])
            sql_string += passed_pose_str
        sql_string = sql_string.rstrip("OR ")
        sql_string += ")"

    if sql_string.endswith("AND "):
        sql_string = sql_string.rstrip("AND ")
    if sql_string.endswith("OR "):
        sql_string = sql_string.rstrip("OR ")
    return sql_string

def write_interaction_index_filtering_str(interaction_list):
    """takes list of interaction info for a given ligand, looks up corresponding interaction index"""
    interaction_info = ["interaction_type", "rec_chain", "rec_resname", "rec_resid", "rec_atom"]
    sql_string = """SELECT interaction_id FROM Interaction_indices WHERE """

    for i in range(4):
        item = interaction_list[i]
        column_name = interaction_info[i]
        if item != "":
            sql_string += "{column} LIKE '%{value}%' AND ".format(column = column_name, value = item)

    sql_string = sql_string.rstrip("AND ")
    return sql_string

def write_interaction_filtering_str(interaction_index_list):
    """takes list of interaction indices and searches for ligand ids which have those interactions"""

    sql_string = """SELECT Pose_id FROM Interaction_bitvectors WHERE """

    for index in interaction_index_list:
        add_str = "Interaction_{index_n} = 1 OR ".format(index_n = index)
        sql_string += add_str

    sql_string = sql_string.rstrip("OR ")

    print(sql_string)
    return sql_string

def select_from_db(conn, filter_sql_string):
    cur = conn.cursor()
    cur.execute(filter_sql_string)

    rows = cur.fetchall()

    if not parsed_opts.no_print:
        for row in rows:
            print(row)

    cur.close()
    conn.close()

    return rows

def write_ligand_filtering_sql(filters):

    sql_ligand_string = "SELECT LigName, best_binding FROM Ligands WHERE "

    ligand_filters = filters["ligand_filters"]
    substruct_flag = ligand_filters['F'][0].upper()
    for kw in ligand_filters.keys():
        fils = ligand_filters[kw]
        if kw == 'N':
            for name in fils:
                name_sql_str = "LigName LIKE '%{value}%' OR ".format(value = name)
                sql_ligand_string += name_sql_str
        if kw == "S":
            for substruct in fils:
                substruct_sql_str = "ligand_smile LIKE '%{value}%' {flag} ".format(value = substruct, flag = substruct_flag)
                sql_ligand_string += substruct_sql_str

    if sql_ligand_string.endswith("AND "):
        sql_ligand_string = sql_ligand_string.rstrip("AND ")
    if sql_ligand_string.endswith("OR "):
        sql_ligand_string = sql_ligand_string.rstrip("OR ")

    return sql_ligand_string

def write_plot_results_sql():
    sql_string = "SELECT LigName, energies_binding, leff, pose_rank FROM Results"

    return sql_string

def write_statevar_sql_call(ligand_name):
    sql_string = "SELECT trans_x, trans_y, trans_z, axisangle_x, axisangle_y, axisangle_z, axisangle_w, dihedrals FROM Results WHERE LigName LIKE '%{value}%' ".format(value = ligand_name)

    return sql_string

def write_input_pdbqt_sql_call(ligand_name):
    sql_string = "SELECT input_pdbqt FROM Ligands WHERE LigName LIKE '%{value}%' ".format(value = ligand_name)

    return sql_string

#####################################################
############ dlg parsing helper functions ###########
#####################################################

def make_list_chunks(lst):
    """break list of dlgs into chuncks of 1000 dlgs"""
    for i in range(0, len(lst), 1000): 
        yield lst[i:i + 1000] 

def dlg_manager(fname):
    """manager for multiprocessing. Parses dlg to dictionary, opens connection to database, writes ligand row"""
    dlg_dict = ad_parser.parse_dlg_gpu(fname)
    dlg_dict = get_best_cluster_poses(dlg_dict)
    resultsAndInteractions = write_results_interaction_rows(dlg_dict)
    results_rows = resultsAndInteractions[0]
    interaction_rows = resultsAndInteractions[1]
    ligand_row = write_ligand_row(dlg_dict)
    return (results_rows, ligand_row, interaction_rows)

def get_best_cluster_poses(ligand_dict):
    """takes input ligand dictionary, reads run pose clusters, adds "cluster_best_run" entry with the top scoring run for each cluster"""
    top_poses = []
    cluster_dict = ligand_dict["clusters"]
    for cluster in cluster_dict:
        top_poses.append(cluster_dict[cluster][0])
    ligand_dict["cluster_top_poses"] = top_poses
    return ligand_dict

###################################################
############ filtering helper functions ###########
###################################################

def make_results_filter_list(filters):
    """takes filters dictionary from option parser. Output list of tuples to be inserted into sql call string"""

    filters_list = []

    #get property filters
    properties_keys = ['eworst',
    'ebest',
    'leworst',
    'lebest']

    property_filters = filters['properties']
    for key in properties_keys:
        if property_filters[key] != None:
            filters_list.append((key, property_filters[key]))

    #get interaction filters
    interaction_keys = ['V',
    'H',
    'HA',
    'HD',
    'R']

    interaction_filters = filters['interactions']
    for key in interaction_filters:
        if interaction_filters[key] != None:
            filters_list.append((key, interaction_filters[key]))

    #get interaction count filters
    interact_count_filters = filters["interactions_count"]
    for count in interact_count_filters:
        filters_list.append(count)

    return filters_list

###################################################
############ output helper functions ##############
###################################################

def write_log(filtered_results, filtered_ligands):
    log_file = parsed_opts.log

    with open(log_file, 'w') as f:
        f.write("Filtered poses:\n")
        f.write("---------------\n")
        for line in filtered_results:
            f.write(" ".join(map(str,line)))
            f.write("\n")
        f.write("\n")
        f.write("***************\n")
        f.write("\n")
        f.write("Filtered Ligands:\n")
        f.write("-----------------\n")
        if filtered_ligands != None:
            for line in filtered_ligands:
                f.write(line + "\n")

    return

"""def write_pdbqts(filtered_ligands, db_file, export_pdbqt_path):
    #Takes list of filtered ligand poses, writes ligand pdbqts with those poses
    ligand_poses = {}

    for pose in filtered_ligands:
        ligand_name = pose[0]

        #get statevars for each pose, add to dictionary with ligand name as key
        statevar_sql_str = write_statevar_sql_call(ligand_name)
        conn = create_connection(db_file)
        statevars = select_from_db(conn, statevar_sql_str)

        statevrs_list = []
        for pose in statevars:
            statevrs_list.append("State: " + " ".join(map(str, pose)).replace("'", "").replace(",", ""))

        if ligand_name in ligand_poses.keys():
            ligand_poses[ligand_name].append(statevrs_list)
        else:
            ligand_poses[ligand_name] = [statevrs_list]

    #get input pdbqt and make temporary input pdbqt file. Then call write_models_from_states
    for ligand in ligand_poses.keys():
        ligand_sv_strings = ligand_poses[ligand]
        conn = create_connection(db_file)
        select_inpt_pdbqt_str = write_input_pdbqt_sql_call(ligand)
        input_pdbqt = select_from_db(conn, select_inpt_pdbqt_str)[0][0]

        tempfile = ligand.rstrip(".pdbqt")+"_temp.pdbqt"
        with open(tempfile, 'w') as f:
            f.write(input_pdbqt)

        print("Making pdbqt for", ligand_name)
        with open("svTemp.txt", 'w') as f:
            for line in statevrs_list:
                f.write(line)
                f.write("\n")

        sp_call = "pythonsh /home/althea/programs/scripts/write_models_from_states.py -v -l " + tempfile + " -s svTemp.txt -o " + export_pdbqt_path+ligand_name
        sp = subprocess.Popen(["/bin/bash", "-i", "-c", sp_call])
        sp.communicate()
        #os.system("pythonsh /home/althea/programs/scripts/write_models_from_states.py -v -l " + tempfile + " -s svTemp.txt -o " + export_pdbqt_path+ligand_name)
        #TODO: change how script is called
        os.system("rm " + tempfile)
        #os.system("rm svTemp.txt")"""

#######################################################################################
################ MAIN #################################################################
#######################################################################################

# create parser
parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
# populate options menu
for group_class in option_groups:
    for group_name, group_info in group_class.items():
        group_desc = group_info['desc']
        group_args = group_info['args']
        # create new group 
        group = parser.add_argument_group(group_name,group_desc)
        # add options to the group
        for name, args in group_args:
            group.add_argument(name, **args)
# parse options
cmdline_opts = sys.argv[1:]
# hack to allow defining options from a file (argparse does not allow to do it in a clean way)
if "--filters_file" in cmdline_opts:
    idx = cmdline_opts.index('--filters_file')
    cmdline_opts.pop(idx)
    ffile = cmdline_opts.pop(idx)
    cmdline_opts += read_filter_file(ffile)
# add a function here to validate the cmdline (repeated options?)
# validate policy of file > cmdline? (or vice versa?)
parsed_opts = parser.parse_args(cmdline_opts)
#print("PARSED OPTIONS", parsed_opts)
file_sources, filters, output, = process_options(parsed_opts)

write_db_flag = False
input_sql = parsed_opts.sql_db
if input_sql == None:
    write_db_flag = True
#############################################
### Find dlgs and write sql database file ###
#############################################

    #open sql db and create tables
    print("Creating database...")
    conn = create_connection(parsed_opts.output_sql)
    create_results_table(conn)
    conn = create_connection(parsed_opts.output_sql)
    create_ligands_table(conn)
    time1 = time.perf_counter()

    #read dlgs in parallel
    manager = ad_parser.DockingResultManager(sources = file_sources, filters=filters, output=output)

    dlg_file_list = manager.files_pool
    dlg_file_list_chunked = list(make_list_chunks(dlg_file_list))

    print(f'Reading dlgs on {mp.cpu_count()} cores')
    all_interactions = []
    all_interactions_split = []
    num_poses = 0
    interaction_chunks = []
    for dlg_chunk in dlg_file_list_chunked:
        results_array = []
        ligands_array = []
        interaction_rows_list = []
        with mp.Pool() as pool:
            parsed_dlgs = pool.map(dlg_manager, dlg_chunk)
            for dlg in parsed_dlgs:
                results_rows = dlg[0]
                ligand_row = dlg[1]
                interaction_rows = dlg[2]
                for pose in results_rows:
                    results_array.append(pose)
                for pose in interaction_rows:
                    interaction_rows_list.append(pose)
                ligands_array.append(ligand_row)

        #find all unique interactions
        all_interactions_data = get_all_interactions(interaction_rows_list, all_interactions, all_interactions_split)
        all_interactions = all_interactions_data[0]
        all_interactions_split = all_interactions_data[1]
        num_poses += len(interaction_rows_list)
        interaction_chunks.append(interaction_rows_list)

        conn = create_connection(parsed_opts.output_sql)
        insert_results(conn, np.array(results_array))
        conn = create_connection(parsed_opts.output_sql)
        insert_ligands(conn, np.array(ligands_array))

    time2 = time.perf_counter()

    #create interaction tables and insert all interactions
    conn = create_connection(parsed_opts.output_sql)
    create_interaction_index_table(conn)
    conn = create_connection(parsed_opts.output_sql)
    create_interaction_bv_table(conn, len(all_interactions))
    conn = create_connection(parsed_opts.output_sql)
    insert_all_interactions(conn, np.array(all_interactions_split))
    #generate interaction bitvector for each pose
    pose_id_list = range(num_poses)
    pose_id_chunks = list(make_list_chunks(pose_id_list))
    for i in range(len(interaction_chunks)):
        with mp.Pool() as pool:
            write_bv_partial = partial(write_interaction_bitvector, all_interactions=all_interactions)
            print(f'Calculating interaction bitvectors on {mp.cpu_count()} cores')
            interaction_bitvectors = pool.starmap(write_bv_partial, zip(pose_id_chunks[i], interaction_chunks[i]))
            conn = create_connection(parsed_opts.output_sql)
            insert_interaction_BVs(conn, interaction_bitvectors, len(all_interactions))

    time3 = time.perf_counter()

    input_sql = parsed_opts.output_sql
time4 = time.perf_counter()

print("Parsing filters...")
#fetch filters for results table, get sql call
filter_list = make_results_filter_list(filters)
sql_filter_str = write_result_filtering_sql(filter_list)
print(sql_filter_str)
conn = create_connection(input_sql)
filtered_results = select_from_db(conn, sql_filter_str)

if filters['filter_ligands_flag']:
    ligand_filters = write_ligand_filtering_sql(filters)
    conn = create_connection(input_sql)
    filtered_ligands = select_from_db(conn, ligand_filters)
else:
    filtered_ligands = None

time5 = time.perf_counter()

if parsed_opts.log != None:
    write_log(filtered_results, filtered_ligands)

time6 = time.perf_counter()

"""if parsed_opts.export_poses_path != None:
    write_pdbqts(filtered_results, input_sql, parsed_opts.export_poses_path)"""

time7 = time.perf_counter()

#make energy and le histograms
if parsed_opts.plot:
    conn = create_connection(input_sql)
    plot_data_str = write_plot_results_sql()
    plot_data = select_from_db(conn, plot_data_str)

    energies = []
    leffs = []
    for row in plot_data:
        if int(row[3]) == 1: #only one point per ligand, for best binding energy
            energies.append(row[1])
            leffs.append(row[2])
    plot_hist(energies, title = "Best pose energies", showSummary=True, xlab=True, bincount=150)
    plot_hist(leffs, title = "Best pose ligand efficiency", showSummary=True, xlab=True, bincount=150)

    colors = []
    for row in plot_data:
        if int(row[3]) == 1:
            ligand_name = row[0]
            if not filters["filter_ligands_flag"]:
                filtered_ligands = []

            if ligand_name in [row[0] for row in filtered_results] and ligand_name in [row[0] for row in filtered_ligands]:
                colors.append("purple")
            elif ligand_name in [row[0] for row in filtered_ligands]:
                colors.append("red")
            elif ligand_name in [row[0] for row in filtered_results]:
                colors.append("blue")
            else:
                colors.append("black")

    plot = plt.scatter(energies, leffs, c=colors)
    if parsed_opts.filters_file != None:
        fig_name = parsed_opts.filters_file.split(".")[0]+".png"
    else:
        fig_name = "out.png"
    plt.savefig(fig_name)

if write_db_flag:
    print("Time to initialize:", time1-time0, "seconds")
    print("Time to read dlgs:", time2-time1, "seconds")
    print("Time to create interaction bitvectors:", time3-time2, "seconds")
    print("Time to write sql database", time4-time3, "seconds")
    print("Total database write time:", time4-time0, "seconds")
print("Filtering time:", time5-time4, "seconds")
if parsed_opts.log != None:
    print("Time to write log:", time6-time5, "seconds")
#if parsed_opts.export_poses_path != None:
    #print("Time to write pdbqts:", time7-time6, "seconds")

