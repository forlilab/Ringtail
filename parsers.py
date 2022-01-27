import os
import gzip
import numpy as np
from glob import glob
import fnmatch

'''class DockingResultManager(object):
    
    def __init__(self, sources = None, filters=None, output=None):
        """ initialize the manager to populate the files list """
        # define file sources to use
        self._process_sources(sources)

    def _process_sources(self, sources):
        """ process the options for input files (parse dictionary) """
        self.files_pool = []
        if sources['file'] is not None:
            # print("DRM> initialized with %d individual files" % len(sources['file']))
            self.files_pool = sources['file']
        # update the files pool with the all the files found in the path
        if sources['file_path'] is not None:
            # print("DRM> scanning path [%s] (recursive: %s, pattern '%s')" % (sources['file_path']['path'],
            #         str(sources['file_path']['recursive']), sources['file_path']['pattern']))
            self.scan_dir(sources['file_path']['path'], 
                          sources['file_path']['pattern'], 
                          sources['file_path']['recursive'])
        # update the files pool with the files specified in the files list
        if sources['file_list'] is not None:
            # print("DRM> searching for files listed in [%s]" % sources['file_list'])
            self.scan_file_list(sources['file_list'])

    def scan_dir(self, path, pattern, recursive=False):
        """ scan for valid output files in a directory 
            the pattern is used to glob files
            optionally, a recursive search is performed
        """
        print("-Scanning directory [%s] for DLG files (pattern:|%s|)" % (path, pattern))
        files = []
        if recursive:
            path = os.path.normpath(path)
            path = os.path.expanduser(path)
            for dirpath, dirnames, filenames in os.walk(path):
                files.extend(os.path.join(dirpath,f) for f in fnmatch.filter(filenames,'*'+pattern))
        else:
            files = glob(os.path.join(path, pattern))
        print("-Found %d files." % len(files))
        self.files_pool.extend(files)'''

def parse_single_dlg(fname, mode='standard'):
    """ parse an ADGPU DLG file uncompressed or gzipped 
    """
    STD_END = 'DOCKED: ENDMDL'
    STD_KW = 'DOCKED: '
    if mode == 'input':
        STD_END = 'FINAL DOCKED STATE:'
        STD_KW = 'INPUT-LIGAND-PDBQT: '

    INPUT_KW = "INPUT LIGAND PDBQT FILE"
    INPUT_END = "FINAL DOCKED STATE"


    # split the first name/extension
    fname_clean = os.path.basename(fname)
    name, ext = os.path.splitext(fname_clean)
    ext = ext[1:].lower()
    if ext == 'gz':
        open_fn = gzip.open
        # split the second name/extension
        name, ext = os.path.splitext(name)
        ext = ext[1:].lower()
    else:
        open_fn = open

    #intialize lists for pose data
    poses = []
    interactions = []
    scores = []
    intermolecular_energy = []
    vdw_hb_desolv = []
    electrostatic = []
    flex_ligand = []
    flexLigand_flexReceptor = []
    internal_energy = []
    torsion = []
    unbound_energy = []
    pose_about = []
    pose_trans = []
    pose_quarternions = []
    pose_dihedrals = []
    clusters = {}
    pose_interact_count = []
    pose_hb_counts = []

    #read poses
    heavy_at_count = 0
    heavy_at_count_complete = False
    #print("Reading", fname)
    with open_fn(fname, 'rb') as fp:
        inside_pose = False
        inside_input = False
        smile_string = ""
        input_pdbqt = []
        num_interact = 0
        for line in fp.readlines():
            line = line.decode("utf-8")
            #store ligand file name
            if line[0:11] == "Ligand file":
                ligname = line.split(":",1)[1].strip()
            #store smile string
            if "REMARK SMILES" in line:
                smile_string = line.split("REMARK SMILES")[-1]
            #store number of runs
            if "Number of runs:" in line:
                nruns = int(line.split()[3])
                cluster_rmsds = list(range(nruns))
                ref_rmsds = list(range(nruns))
            #store input pdbqt lines
            if INPUT_KW in line:
                inside_input = True
            if INPUT_END in line:
                inside_input = False
            if inside_input == True:
                if line.startswith("INPUT-LIGAND-PDBQT") or line.startswith("INPUT-FLEXRES-PDBQT"):
                    if " UNK " in line: #replace ligand atoms ATOM flag with HETATM
                        line = line.replace("ATOM", "HETATM")
                    input_pdbqt.append(' '.join(line.split()[1:]))

            #store poses in each cluster in dictionary as list of ordered runs 
            if "RANKING" in line:
                cluster_num = line.split()[0]
                run = line.split()[2]
                if cluster_num in clusters:
                    clusters[cluster_num].append(int(run))
                else:
                    clusters[cluster_num] = [int(run)]

                cluster_rmsds[int(run)-1] = float(line.split()[4]) #will be stored in order of runs
                ref_rmsds[int(run)-1] = float(line.split()[5])



            #store pose anaylsis
            elif line[0:9] == "ANALYSIS:":
                if inside_pose== False:
                    # first time inside a pose block
                    inside_pose = True
                    interactions.append({})
                    poses.append([])
                # storing interactions
                line = line.split("ANALYSIS:")[1]
                kw, info = line.split(None, 1)
                info = info.replace("{", "")
                info = info.replace("}", "")
                info = info.replace('"', '')
                interactions[-1][kw.lower()] =  [x.strip() for x in info.split(",")]
                if "COUNT" in line:
                    interact_count = int(line.split()[1])
                    interact_list = [''] * interact_count
                    pose_interact_count.append(str(interact_count))
                else:
                    if "TYPE" in line:
                        hb_count = line.count("H")
                        pose_hb_counts.append(hb_count)
                    """kw, info = line.split(None, 1)
                    info = info.replace("{", "")
                    info = info.replace("}", "")
                    info = info.replace(" ", "")
                    info = info.strip("\n")
                    info = info.split(",")
                    for i in range(interact_count):
                        interact_list[i] = interact_list[i] + info[i] + ":"
                    interactions.append(interact_list)"""


            elif STD_END in line:
                inside_pose = False
                heavy_at_count_complete = True
                if mode == 'input':
                    break
            elif (line[0:len(STD_KW)] == STD_KW) and inside_pose:
                # store the pose raw data
                line = line.split(STD_KW)[1]
                poses[-1].append(line)
                if "Estimated Free Energy of Binding" in line:
                    e = float(line.split()[7])
                    scores.append(e)
                if "Final Intermolecular Energy" in line:
                    e = float(line.split()[6])
                    intermolecular_energy.append(e)
                if "vdW + Hbond + desolv Energy" in line:
                    e = float(line.split()[8])
                    vdw_hb_desolv.append(e)
                if "Electrostatic Energy" in line:
                    e = float(line.split()[4])
                    electrostatic.append(e)
                if "Moving Ligand-Fixed Receptor" in line:
                    e = float(line.split()[5])
                    flex_ligand.append(e)
                if "Moving Ligand-Moving Receptor" in line:
                    e = float(line.split()[5])
                    flexLigand_flexReceptor.append(e)
                if "Final Total Internal Energy" in line:
                    e = float(line.split()[7])
                    internal_energy.append(e)
                if "Torsional Free Energy" in line:
                    e = float(line.split()[6])
                    torsion.append(e)
                if "Unbound System's Energy" in line:
                    e = float(line.split()[6])
                    unbound_energy.append(e)
                #store state variables
                if "NEWDPF about" in line:
                    ind_pose_abt = [float(i) for i in line.split()[3:]]
                    pose_about.append(ind_pose_abt)
                if "NEWDPF tran0" in line:
                    translations = [float(i) for i in line.split()[3:]]
                    pose_trans.append(translations)
                if "NEWDPF axisangle0" in line:
                    quaternion = [float(i) for i in line.split()[3:]]
                    pose_quarternions.append(quaternion)
                if "NEWDPF dihe0" in line:
                    dihedrals = [float(i) for i in line.split()[3:]]
                    pose_dihedrals.append(dihedrals)

                # update heavy atom count
                if heavy_at_count_complete:
                    continue
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    # count heavy atoms
                    if not line[-2] == "HD":
                        heavy_at_count += 1
    
    #rewrite interactions format
    """ligand_interaction_keys = ["type",
    "ligname",
    "ligid",
    "chain",
    "residue",
    "resid",
    "recname",
    "recid"]
    reform_interactions = []
    for pose_interactions_dict in interactions:
        num_interactions = pose_interactions_dict["count"][0]
        pose_interact_count.append(num_interactions)
        interaction_strings_list = []
        for key in ligand_interaction_keys:
            interaction_data = pose_interactions_dict[key]
            if type(interaction_data) == list:
                interaction_data_string = ""
                for interaction in interaction_data:
                    interaction_data_string = interaction_data_string + interaction + ", "
            interaction_strings_list.append(interaction_data_string)
            if key == 'type':
                hb_count = interaction_data_string.count("H")
                pose_hb_counts.append(hb_count)
        interactions_string = ""
        for i in range(int(num_interactions)):
            single_interaction_string = ""
            for line in interaction_strings_list:
                line_list = line.split(",")
                single_interaction_string += line_list[i] + ":"
            interactions_string += single_interaction_string.replace(" ", "") + ", "
        reform_interactions.append(interactions_string)"""

    sorted_idx = np.argsort(scores)
    # sort poses, scores, and interactions
    #TODO: make this into a helper function
    poses = [poses[i] for i in sorted_idx]
    scores = [scores[i] for i in sorted_idx]
    interactions = [interactions[i] for i in sorted_idx]
    intermolecular_energy = [intermolecular_energy[i] for i in sorted_idx]
    vdw_hb_desolv = [vdw_hb_desolv[i] for i in sorted_idx]
    electrostatic = [electrostatic[i] for i in sorted_idx]
    flex_ligand = [flex_ligand[i] for i in sorted_idx]
    flexLigand_flexReceptor = [flexLigand_flexReceptor[i] for i in sorted_idx]
    internal_energy = [internal_energy[i] for i in sorted_idx]
    torsion = [torsion[i] for i in sorted_idx]
    unbound_energy = [unbound_energy[i] for i in sorted_idx]
    pose_about = [pose_about[i] for i in sorted_idx]
    pose_trans = [pose_trans[i] for i in sorted_idx]
    pose_quarternions = [pose_quarternions[i] for i in sorted_idx]
    pose_dihedrals = [pose_dihedrals[i] for i in sorted_idx]
    cluster_rmsds = [cluster_rmsds[i] for i in sorted_idx]
    ref_rmsds = [ref_rmsds[i] for i in sorted_idx]
    #reform_interactions = [reform_interactions[i] for i in sorted_idx]
    pose_interact_count = [pose_interact_count[i] for i in sorted_idx]
    pose_hb_counts = [pose_hb_counts[i] for i in sorted_idx]

    if len(poses)==0 or len(scores)==0 or len(interactions)==0 or len(intermolecular_energy)==0 or len(vdw_hb_desolv)==0 or len(electrostatic)==0 or len(internal_energy)==0 or len(torsion)==0 or len(unbound_energy)==0:
        raise ValueError
    # calculate ligand efficiency and deltas from the best pose
    leff = [ x/heavy_at_count for x in scores]
    delta = [x-scores[0] for x in scores ]

    return {'ligname':ligname,
            'ligand_input_pdbqt':input_pdbqt,
            'ligand_smile_string':smile_string,
            'clusters':clusters,
            'cluster_rmsds':cluster_rmsds,
            'ref_rmsds':ref_rmsds,
            'scores':scores, 
            'leff':leff, 
            'delta':delta,             
            'intermolecular_energy':intermolecular_energy,
            'vdw_hb_desolv':vdw_hb_desolv,
            'electrostatics':electrostatic,
            'flex_ligand':flex_ligand,
            'flexLigand_flexReceptor':flexLigand_flexReceptor,
            'internal_energy':internal_energy,
            'torsional_energy':torsion,
            'unbound_energy':unbound_energy,
            'interactions':interactions,
            'num_interactions': pose_interact_count,
            'num_hb':pose_hb_counts,
            #'poses': poses, #store pose pdbqt
            'sorted_runs': [x+1 for x in sorted_idx],
            'pose_about': pose_about,
            'pose_translations':pose_trans,
            'pose_quarternions':pose_quarternions,
            'pose_dihedrals':pose_dihedrals,
            'fname': fname,
            'accepted':list(range(len(scores)))}
