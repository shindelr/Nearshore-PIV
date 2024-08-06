"""
Robin Shindelman
2024-08-06

Interpolation Cross-Validation. This script uses SciKitLearn and NumPy to 
validate an interpolation model defined by a PIV Julia script.
"""

import numpy as np
# from numpy import nan, ravel
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error

from juliacall import Main as jl
from juliacall import Pkg

Pkg.add("Interpolations")
jl.seval("using Interpolations")
jl.seval('include("./julia_interp.jl")')
jl.seval("using .Interp")

jl_itp = jl.seval("Interp.test_interp")

def cross_validate_jl_itp(data, train_i, test_i):
    """
    Call Julia interpolation function in question and perform interpolation
    on a training data set, then compare the results to the known test set which
    is masked.
    """
    training_set = data.copy()
    # Flatten training set and mask the test index but first set aside true val
    y = np.ravel(training_set)[test_i]
    np.ravel(training_set)[test_i] = np.nan

    itpd = np.array(jl_itp(training_set))
    y_hat = np.ravel(itpd)[test_i]

    return mean_squared_error(y_hat, y)


# Set up fo K-Folds validation
kf = KFold(n_splits=9)
mses = []
data = np.genfromtxt("../../tests/juliaOut/datax.csv", delimiter=',')

for train_i, test_i in kf.split(np.ravel(data)):
    print("yess")






