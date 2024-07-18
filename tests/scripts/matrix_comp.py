#!/usr/bin/env python3
"""
Command line script. Takes two csv files as arguments and compares them for 
equality.
"""

import sys
# from pandas import read_csv
import pandas as pd
import numpy as np

def load_data(args):
    """
    Load csv files into pandas dataframes.
    """
    fp1, fp2 = args[1], args[2]
    if 'julia' in fp1:
        NAME_1 = 'julia'
        NAME_2 = 'mlab'
    else:
        NAME_1 = 'mlab'
        NAME_2 = 'julia'

    df1 = pd.read_csv(fp1, sep=',', header=None)
    df2 = pd.read_csv(fp2, sep=',', header=None)

    return df1, df2, (NAME_1, NAME_2)

def compare_data(df1, df2, names: tuple, im):
    """
    Compare dataframes.
    """
    if im:
        df1 = df1.map(parse_im)
        df2 = df2.map(parse_im)

    return df1.compare(df2, result_names=(names[0], names[1]))

def parse_im(val):
    """
    Parse imaginary numbers. Julia puts "im" with some spaces in its complex
    numbers and mlab puts "i". The values must be parsed in order for pandas to
    compare them accurately.
    """
    if 'NaN' in val:
        return np.nan + np.nan*1j
    if 'im' in val:
        return complex(val.replace('im', 'j').replace(' ', ''))
    if 'i' in val:
        return complex(val.replace('i', 'j'))


if __name__ == '__main__':
    print('=====================\nCompare two matrices\n=====================\n')
    if len(sys.argv) < 2:
        print('ERROR: Please enter exactly two filepaths containing CSVs.\n')
        sys.exit(-1)

    # Parse Flags
    im = False
    if len(sys.argv) > 3:
        if sys.argv[3] == '-im':
            im = True

    df1, df2, names = load_data(sys.argv)
    comp_frame = compare_data(df1, df2, names, im)

    if comp_frame.empty:
        print("Matrices are equivalent!\n")
        sys.exit(0)
    else:
        print(f"{len(comp_frame)} different rows detected:\n")
        print(comp_frame)
