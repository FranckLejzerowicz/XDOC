# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os, time
from os.path import isfile, isdir
import pandas as pd

from XDOC.do import DOC_do
from XDOC.boot import DOC_boot
from XDOC.loess import DOC_loess
from XDOC.ci import DOC_ci


def DOC(
        i_otu: str,
        o_outdir: str,
        m_metadata: str = None,
        p_filter: str = None,
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
        p_cores: int = 1,
        verbose: bool = True
):
    """
    A python wrapper of the R wrapper from https://github.com/Russel88/DOC
    to run the whole DOC analysis
    """

    start = time.clock()

    if p_pair and len(p_pair) != 2:
        raise IOError("There should only be two names in pair")

    if not isfile(i_otu):
        raise IOError("No input table found at %s" % i_otu)

    if verbose:
        print('read')
    otu = pd.read_csv(i_otu, header=0, index_col=0, sep='\t')
    if not isdir(o_outdir):
        os.makedirs(o_outdir)

    # if p_filter:
    #     # Filter / Transform OTU-table
    #     otu = do_filter(otu, p_filter)

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

    end = time.clock()
    print('time:', end - start)

    Final = {
        'DO': Dis_Over[2],
        'LME': LME,
        'LOWESS': LOWESS,
        'NEG': NEG,
        'FNS': FNS,
        'BOOT': LOWES,
        'CI': LCIS
    }
    return Final
