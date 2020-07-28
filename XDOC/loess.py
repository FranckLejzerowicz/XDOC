# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from skmisc.loess import loess


def DOC_loess(
        do: list,
        p_pair: str,
        p_span: float,
        p_degree: float,
        p_family: str,
        p_iterations: int,
        p_surface: str):

    # Subset
    OL = do[0]
    DIS = do[1]
    OL_rows, OL_cols = OL.shape
    # Overlap values for loess prediction
    xs = np.linspace(start=0, stop=1, num=1001)

    # Vectorize
    if not p_pair:
        tril = np.tril_indices(OL_rows, k=-1)
        OL_tri = OL.values[tril]
        DIS_tri = DIS.values[tril]
    else:
        OL_tri = OL.values
        DIS_tri = DIS.values

    # # To data frame
    DF_l = pd.DataFrame({'y': DIS_tri, 'x': OL_tri})
    DF_l = DF_l.loc[~DF_l.isna().any(axis=1)]

    # Lowess
    LOW = loess(y=DF_l.y, x=DF_l.x, span=p_span, degree=p_degree,
                family=p_family, iterations=p_iterations, surface=p_surface)
    xs = [round(x, 3) for x in xs if DF_l.x.min() < x < DF_l.x.max()]
    LOW_pred = LOW.predict(newdata=xs)
    LOW_P = pd.DataFrame({"Overlap": xs, "LOWESS": LOW_pred.values})
    LOW_P = LOW_P.loc[~LOW_P.isna().any(axis=1)]
    return LOW_P

