# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np


def DOC_ci(lowp: pd.DataFrame, p_ci: tuple):
    """
    """
    rjsd = lowp.iloc[:, 1:]
    cis = pd.DataFrame(
        np.percentile(
            rjsd, axis=1, q=[x*100 for x in p_ci]),
        index=list(p_ci),
        columns=rjsd.index).T

    LCI = pd.concat([pd.DataFrame({'Overlap': lowp.iloc[:, 0].values}), cis],
                    axis=1, sort=False)
    LCI = LCI.loc[~LCI.isna().any(axis=1)]
    return LCI
