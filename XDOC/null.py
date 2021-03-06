# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from XDOC.otunull import DOC_otunull
from XDOC.doc import DOC


def DOC_null(
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
        p_nulls: int = 1,
        non_zero: bool = True,
        use_mp: bool = False,
        verbose: bool = True):
    nulls = {}
    for i in range(p_nulls):
        # Make NULL
        otu_null = DOC_otunull(otu, non_zero)
        # Run DOC
        doc_null = DOC(
            otu_null,
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
            use_mp,
            verbose
        )
        nulls['Null.%s' % i] = doc_null

    final_nulls = {}
    for null, finals in nulls.items():
        for table, table_pd_ in finals.items():
            table_pd = table_pd_.copy()
            table_pd['Null_dataset'] = null
            final_nulls.setdefault(table, []).append(table_pd)

    return final_nulls
