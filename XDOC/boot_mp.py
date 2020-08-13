# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import random
import pandas as pd
import numpy as np
import multiprocessing as mp
from skmisc.loess import loess
import statsmodels.formula.api as smf


def ma(x, n):
    wins = np.repeat(1 / n, n)
    smoothed_valid_x = np.convolve(wins, x, mode='valid')
    smoothed_same_x = np.convolve(wins, x, mode='same')
    smoothed_x = [x if x in smoothed_valid_x else np.nan for x in smoothed_same_x]
    print()
    return smoothed_x


def init_worker(llb, ol, dis, x, p_p, p_s, m_a, p_sp, p_d, p_f, p_i, p_su):
    global llboot, OL, DIS, xs, p_pair, p_subr, p_mov_avg, p_span, p_degree, p_family, p_iterations, p_surface
    llboot, OL, DIS, xs, p_pair, p_subr, p_span = llb, ol, dis, x, p_p, p_s, p_sp
    p_mov_avg, p_degree, p_family, p_iterations, p_surface = m_a, p_d, p_f, p_i, p_su


def work(item):
    mp_bootstrap(llboot, OL, DIS, xs, p_pair, p_subr, p_mov_avg, p_span, p_degree, p_family, p_iterations, p_surface, item)


def mp_bootstrap(llboot, OL, DIS, xs, p_pair, p_subr, p_mov_avg, p_span, p_degree, p_family, p_iterations, p_surface, item):
    OL_rows, OL_cols = OL.shape
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
    xs = [x for x in xs if DF_l.x.min() < x < DF_l.x.max()]
    LOW_pred = LOW.predict(newdata=xs)
    LOW_P = pd.DataFrame({"rJSD Boot%s" % item: LOW_pred.values})

    # Data frame for lme (slope)
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
    if sum(slope > 0):
        point = slope[slope > 0].index[-1]
        neg_slope = xs[point]
        # Fns
        Fns = sum([1 for x in OL_tri if x > neg_slope]) / len(OL_tri)
    else:
        neg_slope = np.nan
        Fns = 0
    llboot.append([LOW_P, Est, neg_slope, Fns, item])


def get_boot(
        OL: pd.DataFrame,
        DIS: pd.DataFrame,
        xs: np.ndarray,
        p_r: int,
        p_pair: str,
        p_mov_avg: int,
        p_subr: int,
        p_cores: int,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str,
        use_mp: bool):

    if use_mp:
        m = mp.Manager()
        llboot = m.list()
        p_rs = range(p_r)

        if p_cores:
            if p_cores >= p_r:
                nchunks = 1
                cpus = p_r
            else:
                nchunks = int(p_r / p_cores)
                cpus = p_cores
        else:
            cpus = mp.cpu_count()
            if cpus >= p_r:
                nchunks = 1
                cpus = p_r
            else:
                if cpus <= 6:
                    cpus = 6
                else:
                    cpus = 4
                nchunks = int(p_r / cpus)

        print('number of items:', p_r)
        print('number of procs: %s' % cpus)
        print('number of iters:', nchunks)

        p = mp.Pool(initializer=init_worker,
                    initargs=(llboot, OL, DIS, xs, p_pair, p_subr, p_mov_avg, p_span,
                              p_degree, p_family, p_iterations, p_surface),
                    processes=cpus)
        # for _ in tqdm.tqdm(p.imap_unordered(work, p_rs, chunksize=nchunks)):
        #     pass
        for idx, _ in enumerate(p.imap_unordered(work, p_rs, chunksize=nchunks)):
            sys.stdout.write('\rprogress {0:%}'.format(round(idx/p_r, 1)))
        p.close()
        p.join()
        print()

        llboot = [ll[:-1] for ll in sorted(llboot, key=lambda x: x[-1])]

    else:
        llboot = []
        p_rs = range(p_r)

        print("p_rs")
        print(p_rs[:5])
        print("p_rs")
        print(p_rs[-5:])

        for r in p_rs:
            OL_rows, OL_cols = OL.shape
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
            xs = [x for x in xs if DF_l.x.min() < x < DF_l.x.max()]
            LOW_pred = LOW.predict(newdata=xs)
            LOW_P = pd.DataFrame({"rJSD Boot%s" % r: LOW_pred.values})

            # Data frame for lme (slope)
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
            if sum(slope > 0):
                point = slope[slope > 0].index[-1]
                neg_slope = xs[point]
                # Fns
                Fns = sum([1 for x in OL_tri if x > neg_slope]) / len(OL_tri)
            else:
                neg_slope = np.nan
                Fns = 0
            llboot.append([LOW_P, Est, neg_slope, Fns, r])

        llboot = [ll[:-1] for ll in sorted(llboot, key=lambda x: x[-1])]

    print()
    print("llboot[0].iloc[:5, :]")
    print(llboot[0].iloc[:5, :])
    print(llboot[0].iloc[-5:, :])

    print()
    print("llboot[1].iloc[:5, :]")
    print(llboot[1].iloc[:5, :])
    print(llboot[1].iloc[-5:, :])

    print()
    print("llboot[2].iloc[:5, :]")
    print(llboot[2].iloc[:5, :])
    print(llboot[2].iloc[-5:, :])

    print()
    print("llboot[3].iloc[:5, :]")
    print(llboot[3].iloc[:5, :])
    print(llboot[3].iloc[-5:, :])

    return llboot

