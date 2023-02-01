# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import time
import pandas as pd

from XDOC.do import DOC_do
from XDOC.do_mp import DOC_do_mp
from XDOC.boot import DOC_boot
from XDOC.loess import DOC_loess
from XDOC.ci import DOC_ci


def DOC(
        otu: pd.DataFrame,
        p_r: int = 100,
        p_subr: int = 0,
        p_pair: str = None,
        p_mov_avg: int = 5,
        p_ci: tuple = (0.025, 0.5, 0.975,),
        p_span: float = .2,
        p_degree: float = 1.,
        p_family: str = 'symmetric',
        p_iterations: int = 4,
        p_surface: str = 'direct',
        p_cpus: int = 1,
        use_mp: bool = False,
        verbose: bool = True
):
    """
    A python wrapper of the R wrapper from https://github.com/Russel88/DOC
    to run the whole DOC analysis
    """

    # Normalize OTU-table
    otun = otu / otu.sum()

    if verbose:
        print('Dissimilarity and Overlap')
    # Dissimilarity and Overlap
    # start = time.clock()
    # dis_over = DOC_do(otun, p_pair)
    if use_mp:
        dis_over = DOC_do_mp(otun, p_pair, p_cpus)
    else:
        dis_over = DOC_do(otun, p_pair)
    # end = time.clock()
    # print(':' * 30)
    # print('time:', end - start)
    # print(':' * 30)

    if verbose:
        print('Bootstrap lowess and lme')
    # Bootstrap lowess and lme
    BOOT, LME, NEG, FNS = DOC_boot(
        dis_over,
        p_r,
        p_subr,
        p_pair,
        p_mov_avg,
        p_span,
        p_degree,
        p_family,
        p_iterations,
        p_surface,
        p_cpus,
        use_mp
    )

    if verbose:
        print('Lowess confidence intervals')
    # LOWESS CI
    LCIS = DOC_ci(BOOT, p_ci)

    # LOWESS no bootstrap
    if verbose:
        print('Lowess no bootsrap')
    LOWESS = DOC_loess(
        dis_over, p_pair, p_span, p_degree, p_family, p_iterations, p_surface)

    final = {
        'DO': dis_over[2],
        'LME': LME,
        'LOWESS': LOWESS,
        'NEG': NEG,
        'FNS': FNS,
        'BOOT': BOOT,
        'CI': LCIS
    }
    return final
