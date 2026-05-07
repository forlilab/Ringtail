#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ringtail clustering manager
#

import numpy as np
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


class MorganFingerprintCluster:
    """
    Class to handle clustering based on Morgan fingerprints
    Attributes:
        unclustered_data: identifiers for the data being clustered
        rating_data: data used to identify the top scoring/representative item per cluster
        cutoff: cutoff for each cluster
        rdmols: of the data, will be used to create the bitvectors
    """

    def __init__(
        self,
        unclustered_items: iter,
        rating_data: iter,
        rdmols: iter,
        cutoff_distance: float = 0.5,
    ):
        self.unclustered_items = unclustered_items
        self.rating_data = rating_data
        self.rdmols = rdmols
        self.cutoff = cutoff_distance

    def cluster(self) -> tuple[list, list]:
        """
        Clusters self.unclustered_items using morgan finger prints

        Returns:
            list: of clusters
            list: of top rated item for each cluster
        """
        cluster_indices = butina_cluster_fingerprints(
            self.generate_morgan_fingerprints(),
            self.cutoff,
        )
        clusters = []
        for cluster_group in cluster_indices:
            cluster = []
            for index in cluster_group:
                cluster.append(self.unclustered_items[index])
            clusters.append(cluster)

        return clusters, top_score_per_cluster(
            cluster_indices, self.rating_data, self.unclustered_items
        )

    def generate_morgan_fingerprints(
        self,
    ) -> list[DataStructs.cDataStructs.ExplicitBitVect]:
        """
        Generates morgan fingerprints from a serialized rd Mol

        Returns:
            list[DataStructs.cDataStructs.ExplicitBitVect]: list of morgan fingerprints as bitvectors
        """
        mfpgenerator = rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=1024,
        )

        mfps = []
        # prepare each mol json to fingerprint
        for rdmol in self.rdmols:
            # deserialize mols
            mol = Chem.Mol(rdmol)
            Chem.SanitizeMol(mol)
            mfps.append(mfpgenerator.GetFingerprint(mol))

        return mfps


class InteractionBitvectorCluster:
    """
    Class to handle clustering based on ligand-receptor interaction bitvectors
    Attributes:
        unclustered_data: identifiers for the data being clustered
        rating_data: data used to identify the top scoring/representative item per cluster
        cutoff: cutoff for each cluster
        bitvectors: created based on interactions
    """

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
        """
        Clusters self.unclustered_items using interaction bitvectors

        Returns:
            list: of clusters
            list: of top rated item for each cluster
        """
        cluster_indices = butina_cluster_fingerprints(
            [
                DataStructs.CreateFromBitString(bitvector)
                for bitvector in self.bitvectors
            ],
            self.cutoff,
        )
        clusters = []
        for cluster_group in cluster_indices:
            cluster = []
            for index in cluster_group:
                cluster.append(self.unclustered_items[index])
            clusters.append(cluster)

        return clusters, top_score_per_cluster(
            cluster_indices, self.rating_data, self.unclustered_items
        )


def top_score_per_cluster(
    clusters: list[list], rating_data: list, unclustered_items: list
) -> list[int]:
    """
    Finds the stop scored representative for each cluster, based on the supplied
    rating data

    Args:
        clusters (list[list]): list of clusters
        rating_data (list): data used to rate top scoring item
        unclustered_items (list): original data from before clustering

    Returns:
        list[int]: list of items that are representative (best rating) for each cluster
    """

    cluster_representatives = []

    for cluster in clusters:
        cluster_scores = np.array(
            [rating_data[cluster_element] for cluster_element in cluster]
        )
        best_rep_score = unclustered_items[cluster[np.argmin(cluster_scores)]]
        cluster_representatives.append(str(best_rep_score))

    return cluster_representatives


def butina_cluster_fingerprints(fps, cutoff) -> tuple[tuple]:
    """
    Uses Butina algorithm to cluster the provided bitvectors/fingerprints
    https://macinchem.org/2023/03/05/options-for-clustering-large-datasets-of-molecules/

    Args:
        fps (): fingerprints
        cutoff distance (float)
    """
    from rdkit.ML.Cluster import Butina

    # first generate the distance matrix:
    dists = []
    nfps = len(fps)

    n_dists = nfps * (nfps - 1) // 2
    dists = np.empty(n_dists, dtype=np.float64)
    idx = 0
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        n = i
        dists[idx : idx + n] = 1.0 - np.array(sims)
        idx += n

    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs
