# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np


def DOC_ci(BOOT: pd.DataFrame, p_ci: tuple):
    """
    """
    # print()
    # print()
    # print("BOOT")
    # print(BOOT)

    rjsd = BOOT.iloc[:, 1:]
    # print()
    # print()
    # print("rjsd")
    # print(rjsd)

    cis = pd.DataFrame(
        np.percentile(
            rjsd, axis=1, q=[x*100 for x in p_ci]),
        index=list(p_ci),
        columns=rjsd.index).T
    # print()
    # print()
    # print("cis")
    # print(cis)

    # print()
    # print()
    # print("cis.loc[~cis.isna().any(axis=1)]")
    # print(cis.loc[~cis.isna().any(axis=1)])

    LCI = pd.concat([pd.DataFrame({'Overlap': BOOT.iloc[:, 0].values}), cis],
                    axis=1, sort=False)
    # print()
    # print()
    # print("LCI")
    # print(LCI)
    LCI = LCI.loc[~LCI.isna().any(axis=1)]

    # print()
    # print()
    # print("LCI")
    # print(LCI)

    return LCI
