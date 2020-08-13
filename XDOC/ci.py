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

    print("rjsd.loc[:5, :]")
    print(rjsd.loc[:5, :])
    print("rjsd.loc[-5:, :]")
    print(rjsd.loc[-5:, :])

    cis = pd.DataFrame(
        np.percentile(
            rjsd, axis=1, q=[x*100 for x in p_ci]),
        index=list(p_ci),
        columns=rjsd.index).T

    print("cis.loc[:5, :]")
    print(cis.loc[:5, :])
    print(cis.loc[-5:, :])

    LCI = pd.concat([pd.DataFrame({'Overlap': lowp.iloc[:, 0].values}), cis],
                    axis=1, sort=False)
    LCI = LCI.loc[~LCI.isna().any(axis=1)]

    print("LCI.loc[:5, :]")
    print(LCI.loc[:5, :])
    print(LCI.loc[-5:, :])

    return LCI
