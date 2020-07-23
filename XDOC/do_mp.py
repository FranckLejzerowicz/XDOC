# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import tqdm
import pandas as pd
import numpy as np
import itertools

import multiprocessing as mp
from XDOC.rjds import DOC_rjsd


def init_worker(Mat_O, Mat_r, ot):
    global Mat_Overlap_dict, Mat_rJSD_dict, otu
    Mat_Overlap_dict, Mat_rJSD_dict, otu = Mat_O, Mat_r, ot


def work(item):
    mp_do(otu, Mat_Overlap_dict, Mat_rJSD_dict, item)


def mp_do(otu, Mat_Overlap_dict, Mat_rJSD_dict, item):

        i, j = item
        A = otu[[i, j]].copy()
        # Shared species
        shared = (A.astype(bool).sum(axis=1) == 2)
        # Overlap
        x = A.loc[shared, i]
        y = A.loc[shared, j]

        overlap = sum(0.5 * (x + y))

        # Renormalize
        renorm_i = x / sum(x)
        renorm_j = y / sum(y)

        rootJSD = DOC_rjsd(renorm_i, renorm_j)

        # Insert in Matrices
        Mat_Overlap_dict[(i, j)] = overlap
        Mat_rJSD_dict[(i, j)] = rootJSD


def DOC_do_mp(otu: pd.DataFrame, pair: str, p_cores: int):

    cols = otu.columns.tolist()
    samples = len(cols)

    n_pairs = (samples * (samples - 1)) / 2.

    m = mp.Manager()
    Mat_Overlap_d = m.dict()
    Mat_rJSD_d = m.dict()

    # ns = m.Namespace()
    # ns.Mat_Overlap = pd.DataFrame(
    #     [[np.nan] * samples] * samples,
    #     index=cols, columns=cols)
    # ns.Mat_rJSD = pd.DataFrame(
    #     [[np.nan] * samples] * samples,
    #     index=cols, columns=cols)
    print('number of items:', n_pairs)

    if p_cores:
        print('number of procs:', p_cores)
        nchunks = int(n_pairs / p_cores)
    else:
        cpus = mp.cpu_count()
        print('number of procs:', cpus)
        if cpus >= 6:
            nchunks = int(n_pairs / 6)
        else:
            nchunks = int(n_pairs / 2)
    print('number of iters:', nchunks)
    # use all available CPUs
    iter_items = itertools.combinations(cols, 2)
    if p_cores:
        p = mp.Pool(initializer=init_worker, initargs=(Mat_Overlap_d, Mat_rJSD_d, otu), processes=p_cores)
    else:
        p = mp.Pool(initializer=init_worker, initargs=(Mat_Overlap_d, Mat_rJSD_d, otu), processes=1)
    for idx, _ in enumerate(p.imap_unordered(work, iter_items, chunksize=nchunks)):
        sys.stdout.write('\rprogress {0:%}'.format(round(idx/n_pairs, 1)))
    # for _ in tqdm.tqdm(p.imap_unordered(work, iter_items, chunksize=nchunks)):
    #     pass
    p.close()
    p.join()
    print()
    Mat_Overlap = pd.DataFrame(
        [[np.nan] * samples] * samples,
        index=cols, columns=cols)
    Mat_rJSD = pd.DataFrame(
        [[np.nan] * samples] * samples,
        index=cols, columns=cols)
    for (i, j), v in Mat_Overlap_d.items():
        Mat_Overlap.loc[i, j] = v
        Mat_rJSD.loc[i, j] = Mat_rJSD_d[(i, j)]

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

    return List
