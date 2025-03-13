#!/bin/bash

export JULIA_NUM_THREADS=4

# Arguments
N=2
crop_factor="24, 2424, 1, 2048"
final_win_size=16
ol=0.5
out_dir="tests/pipeline_utility_testing/SVSout_23227179_1724441851/pivframes/"

# SINGLE BATCH
in_path="tests/pipeline_utility_testing/testbatches/tmp.eT30jRtvHl.txt"

# MULTI BATCH
# in_dir="tests/pipeline_utility_testing/testbatches/"
# multi_batch=1

verbose=1

# Path to julia script
exec_path="/Users/robinshindelman/repos/Nearshore-Research/juliaPIV/src/PIVPipelineUtility.jl"

# Run julia script
# julia $exec_path $N $crop_factor $final_win_size $ol $out_dir $in_dir $verbose $multi_batch
julia --project=. $exec_path $N $crop_factor $final_win_size $ol $out_dir $in_path $verbose