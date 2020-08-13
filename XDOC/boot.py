# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np
from XDOC.boot_mp import get_boot


def DOC_boot(
        do: list,
        p_r: int,
        p_subr: int,
        p_pair: str,
        p_mov_avg: int,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str,
        p_cpus: int,
        use_mp: bool
):
    """
    """
    # Subset and margin names
    OL = do[0]
    DIS = do[1]

    OL.index = range(1, (OL.shape[0] + 1))
    DIS.index = range(1, (DIS.shape[0] + 1))
    OL.columns = range(1, (OL.shape[1] + 1))
    DIS.columns = range(1, (DIS.shape[1] + 1))

    # Overlap values for loess prediction
    xs = np.array([round(x, 4) for x in np.linspace(start=0, stop=1, num=1001)])

    print("Running bootstraps")
    llboot = get_boot(OL, DIS, xs, p_r, p_pair, p_mov_avg, p_subr, p_cpus,
                      p_span, p_degree, p_family, p_iterations, p_surface, use_mp)

    # Extract and bind lowess, lme, negative slope and Fns separately
    LOWES = pd.concat([x[0] for x in llboot], axis=1, sort=False)
    print()
    print("LOWES")
    print(LOWES.iloc[:5, :])
    print(LOWES.iloc[-5:, :])

    LOWES = pd.concat([pd.DataFrame({'Overlap': [round(float(x), 4) for x in xs]}), LOWES], axis=1, sort=False)

    print()
    print("LOWES")
    print(LOWES.iloc[:5, :])
    print(LOWES.iloc[-5:, :])
    print(LOWESfds)

    LME = pd.DataFrame({'Slope': [x[1] for x in llboot]})
    NEG = pd.DataFrame({'Neg_Slope': [x[2] for x in llboot if str(x[2])!='nan']})
    FNS = pd.DataFrame({'Fns': [x[3] for x in llboot]})

    return [LOWES, LME, NEG, FNS]
