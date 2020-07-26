# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from os.path import isfile, isdir

import pandas as pd
from XDOC.doc import DOC
from XDOC.null import DOC_null
from XDOC.filter import DOC_filter


def xdoc(
        i_otu: str,
        o_outdir: str,
        m_metadata: str = None,
        p_column: str = None,
        p_column_value: tuple = None,
        p_column_quant: int = 0,
        p_filter_prevalence: float = 0,
        p_filter_abundance: float = 0,
        p_filter_order: str = 'meta-filter',
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
        p_nulls: int = 1,
        non_zero: bool = True,
        null: bool = False,
        verbose: bool = True):
    """
    A python wrapper of the R wrapper
    from https://github.com/Russel88/DOC
    to run the whole DOC analysis
    """

    if p_pair and len(p_pair) != 2:
        raise IOError("There should only be two names in pair")

    if not isfile(i_otu):
        raise IOError("No input table found at %s" % i_otu)

    if verbose:
        print('read')
    otu = pd.read_csv(i_otu, header=0, index_col=0, sep='\t')
    if not isdir(o_outdir):
        os.makedirs(o_outdir)

    message = 'input'
    if m_metadata and p_column and p_column_value or p_filter_prevalence or p_filter_abundance:
        # Filter / Transform OTU-table
        otu = DOC_filter(otu, m_metadata, p_filter_prevalence,
                         p_filter_abundance, p_filter_order,
                         p_column, p_column_value, p_column_quant)
        message = 'filtered'

    if otu.shape[0] < 10:
        raise IOError('Too few features in the %s table' % message)

    if verbose:
        print('Table dimension:', otu.shape)

    Final = DOC(
        otu,
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
        p_cpus,
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
        final_nulls = DOC_null(
            otu,
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
            p_cpus,
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
