#!/usr/bin/env python3
"""
Command line script. Takes two csv files as arguments and compares them for 
equality.
"""

import sys
import pandas as pd
import numpy as np

def load_data(filepaths):
    """
    Load csv files into pandas dataframes.
    """
    fp1, fp2 = filepaths[0], filepaths[1]
    if 'julia' in fp1:
        name_1 = 'julia'
        name_2 = 'mlab'
    else:
        name_1 = 'mlab'
        name_2 = 'julia'

    data_frame_1 = pd.read_csv(fp1, sep=',', header=None)
    data_frame_2 = pd.read_csv(fp2, sep=',', header=None)

    return data_frame_1, data_frame_2, (name_1, name_2)

def process_data(data_frame_1, data_frame_2, lang_names: tuple, feat_flags: list[str]):
    """
    Compare dataframes.
    """
    if '-im' in feat_flags:
        data_frame_1 = data_frame_1.map(parse_im)
        data_frame_2 = data_frame_2.map(parse_im)

    if '-r' in feat_flags:
        data_frame_1 = data_frame_1.map(round_im)
        data_frame_2 = data_frame_2.map(round_im)

    if '-s' in feat_flags:
        data_frame_1.to_csv(f"{lang_names[0]}.csv")
        data_frame_2.to_csv(f"{lang_names[1]}.csv")

    return data_frame_1.compare(data_frame_2, result_names=(lang_names[0], lang_names[1]))

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

def round_im(val):
    """
    Round complex numbers to 11 significant digits. Calling parse_im() is a 
    prequisite to this function.
    """
    return round(val.real, 11) + round(val.imag, 11) * 1j


if __name__ == '__main__':
    print('=====================\nCompare two matrices\n=====================\n')
    if len(sys.argv) < 2:
        print('ERROR: Please enter exactly two filepaths containing CSVs.')
        sys.exit(-1)
    if len(sys.argv[1]) < 4 or len(sys.argv[2]) < 4:
        print('Flags must be entered after the filepaths.\n')
        sys.exit(-1)

    # Parse Flags
    # -im       Parses complex numbers
    # -s        Saves dataframes as csv
    # -r        Round to 11 significant figures
    flags = [str(elem) for elem in sys.argv[3:]]

    df1, df2, names = load_data(sys.argv[1:3])
    comp_frame = process_data(df1, df2, names, flags)

    if comp_frame.empty:
        print("Matrices are equivalent!\n")
        sys.exit(0)
    else:
        print(f"{len(comp_frame)} different rows detected:\n")
        print(comp_frame)
