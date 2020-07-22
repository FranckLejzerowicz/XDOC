# ----------------------------------------------------------------------------
# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the MIT License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np


def DOC_rjsd(x, y):
    """
    Jensen-Shannon divergence between two numeric vectors
    """
    z = 0.5 * (x + y)
    rJS = np.sqrt(0.5 * (np.nansum(x * np.log(x / z)) + np.nansum(y * np.log(y / z))))
    return rJS
