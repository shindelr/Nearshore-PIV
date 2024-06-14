#!/usr/bin/env python3
"""
Command line script. Takes two csv files as arguments and compares them for 
equality.
"""

import sys
from pandas import read_csv


print('=====================\nCompare two matrices\n=====================\n')
if len(sys.argv) < 2:
    print('ERROR: Please enter exactly two filepaths containing CSVs.\n')
    sys.exit(-1)

fp1, fp2 = sys.argv[1], sys.argv[2]
if 'julia' in fp1:
    NAME_1 = 'julia'
    NAME_2 = 'mlab'
else:
    NAME_1 = 'mlab'
    NAME_2 = 'julia'

df1 = read_csv(fp1, sep=',', header=None)
df2 = read_csv(fp2, sep=',', header=None)

comp_frame = df1.compare(df2, result_names=(NAME_1, NAME_2))
if comp_frame.empty:
    print("Matrices are equivalent!\n")
    sys.exit(0)
else:
    print(f"{len(comp_frame)} differences detected:\n")
    print(comp_frame)
