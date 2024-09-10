using Base.Threads
using FileIO
using Images
include("main.jl")
using .JuliaPIV

function get_raw_images(PATH::String)::Vector{String}
    files::Vector{String} = readlines(PATH)
    prefix_dir = "tests/pipeline_utility_testing/"
    # Get raw images and prepend the test directory
    return ["$prefix_dir$file" for file in files]
end

function crop_and_pair_images(images::Vector{String}, crop_factor::Tuple{Int, Int, Int, Int})
    # cropped = Vector{Matrix{Gray{N0f8}}}()
    # cropped_pairs = Vector{Tuple{Matrix{Gray{N0f8}}, Matrix{Gray{N0f8}}}}()
    cropped = Vector()
    cropped_pairs = Vector()
    for file_path in images
        img::Matrix{Gray{N0f8}} = Gray.(load(file_path))
        img = img[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]

        @show eltype(img)
        @assert eltype(img) == Matrix{Gray{N0f8}}
        push!(cropped, img)
    end

    display(cropped)

    # for i in 1:length(cropped)-1
    #     push!(cropped_pairs, (cropped[i], cropped[i+1]))
    # end

    return cropped_pairs
end

function io_main(N::T, crop_factor::Tuple{T, T, T, T}, final_win_size::T, 
                ol::Float64, out_dir::String, in_path::String) where {T}
    images = get_raw_images(in_path)
    cropped_pairs = crop_and_pair_images(images, crop_factor)
    @show cropped_pairs[1]
    # @show JuliaPIV.julia_main(cropped_pairs[1])
end

# Testing:
in_args = (
    N = 2,
    crop_factor = (24, 2424, 1, 2048),
    final_win_size = 16,
    ol = 0.5,
    out_dir = "tests/SVSout_23227179_1724441851/pivframes",
    in_dir = "tests/pipeline_utility_testing/testbatches/tmp.2YNbHiPCwK.txt"
)

io_main(in_args[1], in_args[2], in_args[3], in_args[4], in_args[5], in_args[6])