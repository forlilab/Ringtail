import os
import gzip
import numpy as np
from glob import glob
import fnmatch

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
    pose_coordinates = []
    flexible_residues = []
    flexible_res_coords = []

    #read poses
    heavy_at_count = 0
    heavy_at_count_complete = False
    with open_fn(fname, 'rb') as fp:
        inside_pose = False
        inside_input = False
        inside_res = False
        analysis_flag = False
        smile_string = ""
        input_pdbqt = []
        index_map = []
        h_parents = []
        num_interact = 0
        for line in fp.readlines():
            line = line.decode("utf-8")
            #store ligand file name
            if line[0:11] == "Ligand file":
                ligname = line.split(":",1)[1].strip()
            #store smile string
            if "REMARK SMILES" in line and "IDX" not in line:
                smile_string = line.split("REMARK SMILES")[-1]
            #store flexible residue identities
            if "INPUT-FLEXRES-PDBQT: BEGIN_RES" in line:
                flexible_residues.append(line.split()[2])
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
                if line.startswith("INPUT-LIGAND-PDBQT: REMARK SMILES IDX"):
                    index_map += line.lstrip("INPUT-LIGAND-PDBQT: REMARK SMILES IDX" ).rstrip("\n").split()
                if line.startswith("INPUT-LIGAND-PDBQT: REMARK H PARENT"):
                    h_parents += line.lstrip("INPUT-LIGAND-PDBQT: REMARK H PARENT").rstrip("\n").split()

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

            #set make new flexible residue list if we are in the coordinates for a flexible residue
            if "DOCKED: BEGIN_RES" in line:
                flexible_res_coords[-1].append([])
                inside_res = True
            if "DOCKED: END_RES" in line:
                inside_res = False

            #store pose anaylsis
            elif line[0:9] == "ANALYSIS:":
                analysis_flag = True
                if inside_pose== False:
                    # first time inside a pose block
                    inside_pose = True
                    interactions.append({})
                    poses.append([])
                    pose_coordinates.append([])
                    flexible_res_coords.append([])
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

            elif STD_END in line:
                inside_pose = False
                heavy_at_count_complete = True
                if mode == 'input':
                    break
            elif (line[0:len(STD_KW)] == STD_KW) and inside_pose:
                # store the pose raw data
                line = line.split(STD_KW)[1]
                print(line)
                (quit)
                poses[-1].append(line)
                #store pose coordinates
                if "ATOM" in line:
                    if inside_res:
                        line_split = line.split()
                        if "H" not in line_split[2] or "CH" in line_split[2]: #don't save flex res hydrogen coordinates
                            flexible_res_coords[-1][-1].append(line)
                    else:
                        pose_coordinates[-1].append(line.split()[5:8])
                #store pose data
                if "Estimated Free Energy of Binding" in line:
                    try:
                        e = float(line.split()[7])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[6].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    scores.append(e)
                if "Final Intermolecular Energy" in line:
                    try:
                        e = float(line.split()[6])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[5].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    intermolecular_energy.append(e)
                if "vdW + Hbond + desolv Energy" in line:
                    try:
                        e = float(line.split()[8])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[7].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    vdw_hb_desolv.append(e)
                if "Electrostatic Energy" in line:
                    try:
                        e = float(line.split()[4])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[3].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    electrostatic.append(e)
                if "Moving Ligand-Fixed Receptor" in line:
                    try:
                        e = float(line.split()[5])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[4].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    flex_ligand.append(e)
                if "Moving Ligand-Moving Receptor" in line:
                    try:
                        e = float(line.split()[5])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[4].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    flexLigand_flexReceptor.append(e)
                if "Final Total Internal Energy" in line:
                    try:
                        e = float(line.split()[7])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[6].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    internal_energy.append(e)
                if "Torsional Free Energy" in line:
                    try:
                        e = float(line.split()[6])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[5].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
                    torsion.append(e)
                if "Unbound System's Energy" in line:
                    try:
                        e = float(line.split()[6])
                    except ValueError: #catch off-by-one error if number is next to =
                        try:
                            e = float(line.split()[5].lstrip("="))
                        except ValueError:
                            print("ERROR! Cannot parse {this_line} in {this_file}".format(this_line = line, this_file = fname))
                            raise ValueError
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

    sorted_idx = np.argsort(scores)
    # sort poses, scores, and interactions
    #TODO: make this into a helper function
    poses = [poses[i] for i in sorted_idx]
    pose_coordinates = [pose_coordinates[i] for i in sorted_idx]
    flexible_res_coords = [flexible_res_coords[i] for i in sorted_idx]
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
    pose_interact_count = [pose_interact_count[i] for i in sorted_idx]
    pose_hb_counts = [pose_hb_counts[i] for i in sorted_idx]

    if len(poses)==0 or len(scores)==0 or len(interactions)==0 or len(intermolecular_energy)==0 or len(vdw_hb_desolv)==0 or len(electrostatic)==0 or len(internal_energy)==0 or len(torsion)==0 or len(unbound_energy)==0:
        if not analysis_flag:
            raise RuntimeError("No interaction analysis data in {0}. Rerun AD with interaction analysis".format(fname))
        raise ValueError("Incomplete data in " + fname)
    # calculate ligand efficiency and deltas from the best pose
    leff = [ x/heavy_at_count for x in scores]
    delta = [x-scores[0] for x in scores ]

    return {'ligname':ligname,
            'source_file':fname,
            'ligand_input_pdbqt':input_pdbqt,
            'ligand_index_map':index_map,
            'ligand_h_parents':h_parents,
            'pose_coordinates':pose_coordinates,
            'flexible_res_coordinates':flexible_res_coords,
            'flexible_residues':flexible_residues,
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
            'sorted_runs': [x+1 for x in sorted_idx],
            'pose_about': pose_about,
            'pose_translations':pose_trans,
            'pose_quarternions':pose_quarternions,
            'pose_dihedrals':pose_dihedrals,
            'fname': fname}
