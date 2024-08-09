"""
Robin Shindelman
2024-08-06

Interpolation Cross-Validation. This script uses SciKitLearn and NumPy to 
validate an interpolation model defined by a PIV Julia script.
"""

import numpy as np
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

from juliacall import Main as jl
from juliacall import Pkg

# Pkg.add("Interpolations")
jl.seval("using Interpolations")
jl.seval('include("./julia_interp.jl")')
jl.seval("using .Interp")
jl_itp = jl.seval("Interp.regular_interp")
grid_builder = jl.seval("Interp.build_grids_2")

data = np.genfromtxt("../../tests/juliaOut/datax.csv", delimiter=',')
ys, xs, YI, XI = grid_builder(data)

coarse_points = [(y, x, data[j, i])
                 for j, y in enumerate(ys)
                 for i, x in enumerate(xs)
                 ]

train_p, test_p = train_test_split(coarse_points, test_size=0.2, random_state=20)

train_ys = np.unique([p[0] for p in train_p])
train_xs = np.unique([p[1] for p in train_p])
print(train_xs.shape)

# dim1 = round(len(ys) * 0.8) if len(train_p) % (len(ys) * .8) == 0 else len(ys)
# dim2 = round(len(xs) * 0.8) if len(train_p) % (len(xs) * .8) == 0 else len(xs)
train_vals = np.array([p[2] for p in train_p]).reshape((len(ys), len(xs)))

# fine_grid_vals = jl_itp(train_vals, train_xs, train_ys, XI, YI)
