#!/bin/bash

export JULIA_NUM_THREADS=4

# Arguments
N=3
crop_factor="24, 2424, 1, 2048"
final_win_size=16
ol=0.5
out_dir="tests/pipeline_utility_testing/SVSout_23227179_1724441851/pivframes/"
# in_path="tests/pipeline_utility_testing/testbatches/"
in_path="tests/pipeline_utility_testing/testbatches/tmp.7EOJPExxgd.txt"
# multi_batch=1

# Path to julia script
exec_path="/Users/robinshindelman/repos/Nearshore-Research/juliaPIV/src/io.jl"

# Run julia script
# julia $exec_path $N $crop_factor $final_win_size $ol $out_dir $in_path $multi_batch
julia $exec_path $N $crop_factor $final_win_size $ol $out_dir $in_path