# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def DOC_filter(otu: pd.DataFrame, p_filter_prevalence: float,
               p_filter_abundance: float) -> pd.DataFrame:

    preval, abund = p_filter_prevalence, p_filter_abundance
    if not preval + abund:
        return otu
    otu_filt = otu.copy()
    # get the min number of samples based on prevalence percent
    if preval < 1:
        n_percent = otu.shape[1] * preval
    else:
        n_percent = preval
    # abundance filter in terms of min reads counts
    otu_percent = otu_filt.copy()
    otu_percent_sum = otu_filt.sum(1)
    if abund < 1:
        otu_percent = otu_filt / otu_filt.sum()
        otu_percent_sum = otu_filt.sum(1) / otu_filt.sum(1).sum()

    abund_mode = 'sample'

    # remove features from feature table that are not present
    # in enough samples with the minimum number/percent of reads in these samples
    if abund_mode == 'sample':
        otu_filt = otu_filt.loc[(otu_percent > abund).sum(1) > n_percent, :]
    elif abund_mode == 'dataset':
        otu_filt = otu_filt.loc[otu_percent_sum > abund, :]
    elif abund_mode == 'both':
        otu_filt = otu_filt.loc[(otu_percent > abund).sum(1) > n_percent, :]
        fil_pd_percent_sum = otu_filt.sum(1)
        if abund < 1:
            fil_pd_percent_sum = otu_filt.sum(1) / otu_filt.sum(1).sum()
        otu_filt = otu_filt.loc[fil_pd_percent_sum > abund, :]
    else:
        raise Exception('"%s" mode not recognized' % abund_mode)
    otu_filt = otu_filt.loc[otu_filt.sum(1) > 0, otu_filt.sum(0) > 0]
    return otu_filt
