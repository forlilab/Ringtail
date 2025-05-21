import numpy as np
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


class MorganFingerprintCluster:
    def __init__(
        self,
        unclustered_items: iter,
        rating_data: iter,
        rdmols: iter,
        cutoff_distance: float = 0.5,
    ):
        self.unclustered_items = unclustered_items
        self.rating_data = rating_data
        self.cutoff = cutoff_distance
        self.rdmols = rdmols

    def cluster(self):
        bv_clusters = cluster_fingerprints(
            self.generate_morgan_fingerprints(),
            self.cutoff,
        )

        return bv_clusters, top_leff_per_cluster(
            bv_clusters, self.rating_data, self.unclustered_items
        )

    def generate_morgan_fingerprints(self):
        mfpgenerator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        mfps = []
        # prepare each mol json to fingerprint
        for rdmol in self.rdmols:
            # deserialize mols (JSONToMols returns tuple, hence the [0] subscript)])
            mol = Chem.Mol(rdmol)
            # need to sanitize each mol like this (as opposed to inline) as the method has a return value of 'santize_flag'
            Chem.SanitizeMol(mol)
            mfps.append(mfpgenerator.GetFingerprint(mol))

        return mfps


class InteractionBitvectorCluster:

    def __init__(
        self,
        unclustered_items: iter,
        rating_data: iter,
        bitvectors: dict,
        cutoff_distance: float = 0.5,
    ):
        self.unclustered_items = unclustered_items
        self.rating_data = rating_data
        self.cutoff = cutoff_distance
        self.bitvectors = bitvectors

    def cluster(self):
        bv_clusters = cluster_fingerprints(
            [
                DataStructs.CreateFromBitString(bitvector)
                for bitvector in self.bitvectors
            ],
            self.cutoff,
        )
        return bv_clusters, top_leff_per_cluster(
            bv_clusters, self.rating_data, self.unclustered_items
        )


def top_leff_per_cluster(bv_clusters, rating_data, unclustered_items):
    # interaction clusters representative pose ids
    int_rep_poseids = []

    for cluster in bv_clusters:
        # element 1 in individual pose id item is the ligand efficiency (leff)
        c_leffs = np.array(
            [rating_data[cluster_element] for cluster_element in cluster]
        )
        best_lig_c = unclustered_items[cluster[np.argmin(c_leffs)]]
        int_rep_poseids.append(str(best_lig_c))

    # prepare packet to write to db
    return int_rep_poseids


def cluster_fingerprints(fps, cutoff):
    """
    https://macinchem.org/2023/03/05/options-for-clustering-large-datasets-of-molecules/

    Args:
        fps (): fingerprints
        cutoff distance (float)
    """
    from rdkit.SimDivFilters import rdSimDivPickers
    from rdkit.ML.Cluster import Butina

    lp = rdSimDivPickers.LeaderPicker()
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)

    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    picks = lp.LazyBitVectorPick(fps, nfps, 0.5)
    pickfps = [fps[x] for x in picks]
    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs
