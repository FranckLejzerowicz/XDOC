# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np


def DOC_otunull(otu, non_zero=True):
    if non_zero:
        # Shuffle all non-zero values
        otu_null = []
        for row in otu.values:
            if sum(row == 0):
                idx = np.nonzero(row)
                row[idx] = np.random.permutation(row[idx])
            otu_null.append(row)
        otu_null = pd.DataFrame(otu_null)
    else:
        # Shuffle all values
        otu_null = otu.sample(frac=1, axis=1).sample(frac=1)

    otu_null.columns = otu.columns
    otu_null.index = otu.index
    return otu_null
