# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import abspath, dirname, isfile, isdir, splitext
import pandas as pd
import numpy as np
from scipy.signal import lfilter
import itertools
import random
from skmisc.loess import loess
import statsmodels.formula.api as smf

import altair as alt


def DOC_loess(
        do: list,
        p_pair: str,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str):

    # Subset
    OL = do[0]
    DIS = do[1]
    OL_rows, OL_cols = OL.shape
    # Overlap values for loess prediction
    xs = np.linspace(start=0, stop=1, num=1001)

    # Vectorize
    if not p_pair:
        tril = np.tril_indices(OL_rows, k=-1)
        OL_tri = OL.values[tril]
        DIS_tri = DIS.values[tril]
    else:
        OL_tri = OL.values
        DIS_tri = DIS.values

    # # To data frame
    DF_l = pd.DataFrame({'y': DIS_tri, 'x': OL_tri})
    DF_l = DF_l.loc[~DF_l.isna().any(axis=1)]

    # Lowess
    LOW = loess(y=DF_l.y, x=DF_l.x, span=p_span, degree=p_degree,
                family=p_family, iterations=p_iterations, surface=p_surface)
    LOW_pred = LOW.predict(newdata=xs)
    LOW_P = pd.DataFrame({"Overlap": xs, "LOWESS": LOW_pred.values})
    LOW_P = LOW_P.loc[~LOW_P.isna().any(axis=1)]
    return LOW_P


def DOC_ci(lowp: pd.DataFrame, p_ci: tuple):
    """
    """
    rjsd = lowp.iloc[:,1:]
    cis = pd.DataFrame(
        np.percentile(
            rjsd, axis=1, q=[x*100 for x in p_ci]),
        index=list(p_ci),
        columns=rjsd.index).T

    LCI = pd.concat([pd.DataFrame({'Overlap': lowp.iloc[:, 0].values}), cis],
                    axis=1, sort=False)
    LCI = LCI.loc[~LCI.isna().any(axis=1)]
    return LCI


def ma(x, n):
    wins = np.repeat(1 / n, n)
    smoothed_valid_x = np.convolve(wins, x, mode='valid')
    smoothed_same_x = np.convolve(wins, x, mode='same')
    smoothed_x = [x if x in smoothed_valid_x else np.nan for x in smoothed_same_x]
    return smoothed_x


def get_boot(
        OL: pd.DataFrame,
        DIS: pd.DataFrame,
        xs: np.ndarray,
        p_r: int,
        p_pair: str,
        p_mov_avg: int,
        p_subr: int,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str):
    llboot = []
    OL_rows, OL_cols = OL.shape
    for i in range(p_r):
        if not p_pair:
            # Sample subjects
            if p_subr:
                Samp = random.choices(OL.index, k=p_subr)
            else:
                Samp = random.choices(OL.index, k=OL_rows)
            # Subset
            OL_sub = OL.loc[Samp, Samp]
            DIS_sub = DIS.loc[Samp, Samp]
            # Vectorize
            tril = np.tril_indices(len(Samp), k=-1)
            OL_tri = OL_sub.values[tril]
            DIS_tri = DIS_sub.values[tril]
        else:
            # Sample subjects
            if p_subr:
                Sampr = random.choices(OL.index, k=p_subr)
                Sampc = random.choices(OL.columns, k=p_subr)
            else:
                Sampr = random.choices(OL.index, k=OL_rows)
                Sampc = random.choices(OL.columns, k=OL_cols)
            # Subset
            OL_sub = OL.loc[Sampr, Sampc]
            DIS_sub = DIS.loc[Sampr, Sampc]
            # Vectorize
            OL_tri = OL_sub.T.stack(dropna=False).values
            DIS_tri = DIS_sub.T.stack(dropna=False).values

        # To data frame
        DF_l = pd.DataFrame({'y': DIS_tri, 'x': OL_tri})
        DF_l = DF_l.loc[~DF_l.isna().any(axis=1)]
        # Lowess
        LOW = loess(y=DF_l.y, x=DF_l.x, span=p_span, degree=p_degree,
            family=p_family, iterations=p_iterations, surface=p_surface)
        LOW_pred = LOW.predict(newdata=xs)
        LOW_P = pd.DataFrame({"rJSD Boot%s" % i: LOW_pred.values})
        #
        # # Data frame for lme (slope)
        tril = np.tril_indices(OL_sub.shape[1], k=-1)
        OL_vals = OL_sub.values[tril]
        DIS_vals = DIS_sub.values[tril]
        Tris = pd.DataFrame({
            'Row': OL_sub.columns[tril[1]],
            'Col': OL_sub.index[tril[0]],
            'OL': OL_vals,
            'DIS': DIS_vals
        })
        # Remove data with Overlap below median
        Tris_sub = Tris.loc[Tris['OL'] >= np.nanmedian(OL_sub), :]

        # LME
        md = smf.mixedlm("DIS ~ OL", Tris_sub, groups=Tris_sub["Row"])
        fit = md.fit()
        Est = fit.params['OL']

        ## Detect negative slope
        # Smooth prediction
        low_ma = ma(LOW_pred.values, p_mov_avg)
        slope = pd.Series(low_ma).diff() / pd.Series(xs).diff()
        point = slope[slope > 0].index[-1]
        neg_slope = xs[point-1]

        # Fns
        Fns = sum([1 for x in OL_tri if x > neg_slope]) / len(OL_tri)

        llboot.append([LOW_P, Est, neg_slope, Fns])

    return llboot


def DOC_boot(
        do: list,
        p_r: int,
        p_subr: int,
        p_pair: str,
        p_mov_avg: int,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str,
        p_cores: int
):
    """
    """
    # Subset and margin names
    OL = do[0]
    DIS = do[1]
    
    OL.index = range(1, (OL.shape[0]+1))
    DIS.index = range(1, (DIS.shape[0]+1))
    OL.columns = range(1, (OL.shape[1]+1))
    DIS.columns = range(1, (DIS.shape[1]+1))
    
    # Overlap values for loess prediction
    xs = np.linspace(start=0, stop=1, num=1001)

    # Start parallel
    if p_cores == 1:
        pass
    else:
        pass

    print("Running bootstraps")
    
    # Progress bar
#    pb = txtProgressBar(max = R, style = 3)
#    progress = function(n) setTxtProgressBar(pb, n)
#    opts = list(progress = progress)
    llboot = get_boot(OL, DIS, xs, p_r, p_pair, p_mov_avg, p_subr,
                      p_span, p_degree, p_family, p_iterations, p_surface)
    # llboot.append([LOW_P, Est, neg_slope, Fns])

    # Extract and bind lowess, lme, negative slope and Fns seperately
    LOWES = pd.concat([x[0] for x in llboot], axis=1, sort=False)
    LOWES = pd.concat([pd.DataFrame({'Overlap': xs}), LOWES], axis=1, sort=False)
    LME = pd.DataFrame({'Slope': [x[1] for x in llboot]})
    NEG = pd.DataFrame({'Neg_Slope': [x[2] for x in llboot]})
    FNS = pd.DataFrame({'Fns': [x[3] for x in llboot]})
    return [LOWES, LME, NEG, FNS]
    

def DOC_rjsd(x, y):
    """
    Jensen-Shannon divergence between two numeric vectors
    """
    z = 0.5 * (x + y)
    rJS = np.sqrt(0.5 * (np.nansum(x * np.log(x / z)) + np.nansum(y * np.log(y / z))))
    return rJS


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

        overlap = sum(0.5 * (x + y))

        # Renormalize
        renorm_i = x / sum(x)
        renorm_j = y / sum(y)

        # rJSD
        rootJSD = DOC_rjsd(renorm_i, renorm_j)

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

    return List


def xdoc(
        i_otu: str,
        m_metadata: str,
        o_outdir: str,
        p_r: int,
        p_subr: int,
        p_pair: str,
        p_mov_avg: int,
        p_ci: tuple,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str,
        p_cores: int,
        verbose: bool
):
    """
    A python wrapper of the R wrapper from https://github.com/Russel88/DOC
    to run the whole DOC analysis
    """

    if p_pair and len(p_pair) != 2:
        raise IOError("There should only be two names in pair")

    if not isfile(i_otu):
        raise IOError("No input table found at %s" % i_otu)

    if verbose:
        print('read')
    otu = pd.read_csv(i_otu, header=0, index_col=0, sep='\t')
    # Normalize OTU-table
    otun = otu / otu.sum()

    if verbose:
        print('Dissimilarity and Overlap')
    # Dissimilarity and Overlap
    Dis_Over = DOC_do(otun, p_pair)

    if verbose:
        print('Bootstrap lowess and lme')
    # Bootstrap lowess and lme
    LOWES, LME, NEG, FNS = DOC_boot(
        Dis_Over,
        p_r,
        p_subr,
        p_pair,
        p_mov_avg,
        p_span,
        p_degree,
        p_family,
        p_iterations,
        p_surface,
        p_cores
    )

    if verbose:
        print('Lowess confidence intervals')
    # LOWESS CI
    LCIS = DOC_ci(LOWES, p_ci)

    # LOWESS no bootstrap
    if verbose:
        print('Lowess no bootsrap')
    LOWESS = DOC_loess(Dis_Over, p_pair, p_span, p_degree, p_family, p_iterations, p_surface)

    Final = {
        'DO': Dis_Over[2],
        'LME': LME,
        'LOWESS': LOWESS,
        'NEG': NEG,
        'FNS': FNS,
        'BOOT': LOWES,
        'CI': LCIS
    }
    if verbose:
        print('Writing:')
    if not isdir(o_outdir):
        os.makedirs(o_outdir)
    for table, table_pd in Final.items():
        path = '%s/%s.tsv' % (o_outdir, table)
        if verbose:
            print(' -', path)
        if path == 'DO':
            table_pd.to_csv(path, index=True, sep='\t')
        else:
            table_pd.to_csv(path, index=False, sep='\t')
