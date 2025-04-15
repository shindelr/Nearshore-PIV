#!/bin/python3
"""
A script to run the Julia PIV portion of the pipeline.
Uses 4 threads by default to make the comparison of frame pairs run a little 
quicker.
Meant to be executed from the top level of the juliaPIV directory.
"""

import os
import argparse
import subprocess
import logging
from multiprocessing import Pool

logging.basicConfig(level=logging.INFO, format="%(asctime)s -- %(levelname)s -- %(message)s")
NPROC = 4

def launch_batch(file, args):
    args.in_path = file
    logging.info(f"Processing {args.in_path}")
    run_pipe(args=args)


def batches(abs_batch_dir):
    return [os.path.join(abs_batch_dir, f) for f in os.listdir(abs_batch_dir)]


def run_pipe(args):
    # Nifty way of getting absolute path of the script to run, regardless of cwd
    exec_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'piv_build/bin/PIVPipelineUtility'))
    cmmd = [exec_path,
            str(args.N),
            str(args.crop_factor), 
            str(args.final_win_size), 
            str(args.ol), 
            args.out, 
            args.in_path, 
            str(args.verbosity)]
    subprocess.run(cmmd)


def main():
    output_mat_structure = """
This script will output the results of JuliaPIV as a .mat file with the following
structure:
    
x: [255x299 double]
y: [255x299 double]
pass_sizes: [3x2 double]
    overlap: 0.5
    method: 'multin'
        fn: {list of jpg files}
            u: [255x299 double]
            v: [255x299 double]
        npts: [255x299 double]  # number of data points that weren't NaN prior to time-average
        uStd: [255x299 double]  # standard deviation of the N results
        vStd: [255x299 double]  # ditto
"""

    parser = argparse.ArgumentParser(prog='PIV Pipeline Utility',
                                     description='Run Julia PIV on a batch of frames.',
                                     epilog=output_mat_structure)
    
    parser.add_argument('-N', 
                        help="The number of frames to average together at once.",
                        type=int,
                        default=2)
    
    parser.add_argument('--crop_factor',
                        help='Gives a box to extract from the raw image. Should be a tuple of 4 ints',
                        default="24, 2425, 1, 2048")
    
    parser.add_argument('--final_win_size',
                        help='Final window size to evaluate PIV at.',
                        type=int,
                        default=16)
    
    parser.add_argument('--ol', 
                        help='Window overlap for frame comparison.',
                        type=float,
                        default=0.5)
    
    parser.add_argument('--out',
                        help='Where to output .mat files')
    
    parser.add_argument('--in_path',
                        help='.txt file containing image paths.')
    
    parser.add_argument('--verbosity', 
                        help='1 for verbose print statements, 0 otherwise.',
                        type=int,
                        default=0)

    args = parser.parse_args()

    txt_list = batches(args.in_path)
    logging.info(f"Found {len(txt_list)} .txt files\n")

    with Pool(processes=NPROC) as pool:
        pool.starmap(launch_batch, [(file, args) for file in txt_list])


if __name__ == '__main__':
    main()
