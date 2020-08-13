# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np
import itertools

from scipy.spatial.distance import jensenshannon


def DOC_do(otu: pd.DataFrame, pair: str):
    """
    :param otun:
    :param pair:
    :return: a matrix of Overlap, a matrix of rJSD (Dissimilarities), and a dataframe with Overlap and rJSD as vectors
    """
    cols = otu.columns.tolist()
    samples = len(cols)

    Mat_Overlap = pd.DataFrame(
        [[np.nan] * samples] * samples,
        index=cols, columns=cols
    )
    Mat_rJSD = Mat_Overlap.copy()

    C = (samples * (samples + 1)) / 2.

    for c, (i, j) in enumerate(itertools.combinations(cols, 2)):
        if c % 10000 == 0:
            print('%s / %s' % (c, C))
        A = otu[[i, j]].copy()
        # Shared species
        shared = (A.astype(bool).sum(axis=1) == 2)
        # Overlap
        x = A.loc[shared, i]
        y = A.loc[shared, j]

        overlap = round(sum(0.5 * (x + y)), 5)

        # Renormalize
        renorm_i = x / sum(x)
        renorm_j = y / sum(y)

        # rJSD
        rootJSD = jensenshannon(renorm_i, renorm_j)

        # Insert in Matrices
        Mat_Overlap.loc[i, j] = overlap
        Mat_rJSD.loc[i, j] = rootJSD

    if pair:
        pairv = [pair[0] if pair[0] in x else x for x in cols]
        pairv = [pair[1] if pair[1] in x else x for x in pairv]
        if len(set(pairv)) != 2:
            raise IOError("Names of pairs do not match column names")
        Mat_Overlap.index = pairv
        Mat_Overlap.columns = pairv
        Mat_Overlap = Mat_Overlap.loc[pair[0], pair[1]]

        Mat_rJSD.index = pairv
        Mat_rJSD.columns = pairv
        Mat_rJSD = Mat_rJSD.loc[pair[0], pair[1]]

        DF = pd.concat(
            {'Overlap': Mat_Overlap.T.stack(dropna=False),
             'rJSD': Mat_rJSD.T.stack(dropna=False)}, axis=1
        )
        DF = DF.loc[~DF.Overlap.isna()]

        List = [Mat_Overlap, Mat_rJSD, DF]

    else:

        DF = pd.concat(
            {'Overlap': Mat_Overlap.T.stack(dropna=False),
             'rJSD': Mat_rJSD.T.stack(dropna=False)}, axis=1
        )
        DF = DF.loc[~DF.Overlap.isna()]
        Mat_Overlap_t = Mat_Overlap.T
        Mat_rJSD_t = Mat_rJSD.T

        Mat_Overlap[Mat_Overlap.isna()] = 0
        Mat_Overlap_t[Mat_Overlap_t.isna()] = 0
        Mat_rJSD[Mat_rJSD.isna()] = 0
        Mat_rJSD_t[Mat_rJSD_t.isna()] = 0

        Mat_Overlap_new = Mat_Overlap + Mat_Overlap_t
        Mat_rJSD_new = Mat_rJSD + Mat_rJSD_t

        for col in cols:
            Mat_Overlap_new.loc[col, col] = np.nan
            Mat_rJSD_new.loc[col, col] = np.nan

        List = [Mat_Overlap_new, Mat_rJSD_new, DF]

    print('Mat_Overlap_new')
    print(Mat_Overlap_new.iloc[:5, :5])
    print("Mat_rJSD_new")
    print(Mat_rJSD_new.iloc[:5, :5])
    print("DF")
    print(DF.iloc[:5, :])

    return List
