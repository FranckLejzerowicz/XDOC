# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from XDOC.null import xdoc_null
from XDOC.doc import DOC
# def do_filter(otu: pd.DataFrame, p_filter: str):


def xdoc(
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
        p_nulls: int = 1,
        non_zero: bool = True,
        null: bool = False,
        verbose: bool = True):
    """
    A python wrapper of the R wrapper
    from https://github.com/Russel88/DOC
    to run the whole DOC analysis
    """
    Final = DOC(
        i_otu,
        o_outdir,
        m_metadata,
        p_filter,
        p_r,
        p_subr,
        p_pair,
        p_mov_avg,
        p_ci,
        p_span,
        p_degree,
        p_family,
        p_iterations,
        p_surface,
        p_cores,
        verbose
    )
    if verbose:
        print('Writing:')
    for table, table_pd in Final.items():
        path = '%s/%s.tsv' % (o_outdir, table)
        if verbose:
            print(' -', path)
        if table == 'DO':
            table_pd.to_csv(path, index=True, sep='\t')
        else:
            table_pd.to_csv(path, index=False, sep='\t')

    if null:
        final_nulls = xdoc_null(
            i_otu,
            o_outdir,
            m_metadata,
            p_filter,
            p_r,
            p_subr,
            p_pair,
            p_mov_avg,
            p_ci,
            p_span,
            p_degree,
            p_family,
            p_iterations,
            p_surface,
            p_cores,
            p_nulls,
            non_zero,
            verbose
        )
        for table, table_pds in final_nulls.items():
            table_pd = pd.concat(table_pds)
            path = '%s/null_%s.tsv' % (o_outdir, table)
            if verbose:
                print(' -', path)
            if table == 'DO':
                table_pd.to_csv(path, index=True, sep='\t')
            else:
                table_pd.to_csv(path, index=False, sep='\t')
