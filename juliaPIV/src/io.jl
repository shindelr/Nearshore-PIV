module PIVPipelineUtility
export julia_main 

using Base.Threads
using FileIO
using Images
using Plots
using Statistics
using MAT
using DelimitedFiles
include("main.jl")
using .JuliaPIV

function get_raw_images(PATH::String)::Vector{String}
    files::Vector{String} = readlines(PATH)
    prefix_dir = "tests/pipeline_utility_testing/"
    # Get raw images and prepend the test directory
    return ["$prefix_dir$file" for file in files]
end

function crop_and_pair_images(images::Vector{String}, crop_factor::NTuple{4, Int32})
    cropped = Vector()
    cropped_pairs = Vector()
    for file_path in images
        img::Matrix{Gray{N0f8}} = load(file_path)
        img = img[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        push!(cropped, img)
    end

    for i in 1:length(cropped)-1
        push!(cropped_pairs, (cropped[i], cropped[i+1]))
    end

    return cropped_pairs
end

function statistics_of_piv(piv_results, n::Int32)
    println("Calculating statistics...")
    # Preallocations
    xs = Vector{Matrix{Float32}}()
    ys = Vector{Matrix{Float32}}()
    us = Vector{Matrix{Float32}}()
    vs = Vector{Matrix{Float32}}()
    u_avs = Vector{Matrix{Float32}}()
    u_stds = Vector{Matrix{Float32}}()
    v_avs = Vector{Matrix{Float32}}()
    v_stds = Vector{Matrix{Float32}}()
    x_avs = Vector{Matrix{Float32}}()
    y_avs = Vector{Matrix{Float32}}()
    npts = Vector{Matrix{Int32}}()

    # Unpack results
    for result in piv_results
        push!(xs, result[1][1])
        push!(ys, result[1][2])
        push!(us, result[2][1])
        push!(vs, result[2][2])
    end

    # Calculate statistics
    i = 1
    while i <= length(us)
        # Make groups of N results for averages, stds, etc.
        un_group = us[i:min(i+n-1, length(us))]
        vn_group = vs[i:min(i+n-1, length(us))]
        xn_group = xs[i:min(i+n-1, length(us))]
        yn_group = ys[i:min(i+n-1, length(us))]

        # Create binary masks for each matrix in u-group
        nan_binary_masks = Vector{Matrix{Bool}}([isnan.(u) for u in un_group])
        push!(npts, sum(nan_binary_masks))

        # Compute averages and stds
        push!(u_avs, mean(un_group))
        push!(u_stds, std(un_group))
        push!(v_avs, mean(vn_group))
        push!(v_stds, std(vn_group))
        push!(x_avs, mean(xn_group))
        push!(y_avs, mean(yn_group))

        # Increment by N to get the appropriate num of averages etc.
        i += n
    end
    
    return ((x_avs, y_avs), (u_avs, v_avs), (u_stds, v_stds), npts)
end

function heatmap_mult_images(images::Vector{Matrix{Float32}}, title::String)
    for im in images
        display(plot(heatmap(im, 
                title = title, 
                aspect_ratio = :equal, 
                limits=(0, 200), 
                xlimits=(0, 385))))
    end
end

function parse_image_names(images::Vector{String}, N::Int32)::Vector{String}
    i = 1
    image_names = Vector{String}()
    while i < length(images)
        push!(image_names, images[i][69:end-4])
        i += N
    end
    return image_names
end

# SINGLE BATCH
function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String) where {T}

    # Image pre-processing
    images = get_raw_images(in_path)
    if length(images) ÷ N <= 1
        error("\n\nNumber of images in directory ($(length(images))) not divisible by N ($N).\n\n")
    end
    cropped_pairs = crop_and_pair_images(images, crop_factor)

    # PIV, could thread here?
    raw_piv_results = Vector{Tuple{Tuple{Matrix{T}, Matrix{T}}, Tuple{Matrix{T}, 
                                Matrix{T}}, Vector{Int32}} where {T}}()

    Threads.@threads for pair in cropped_pairs
        push!(raw_piv_results, JuliaPIV.main(pair, Int32(final_win_size), Float32(ol)))
    end

    # npts = count_nan(raw_piv_results[])
    ((x_avs, y_avs), (u_avs, v_avs), (u_stds, v_stds), npts) = statistics_of_piv(raw_piv_results, N)

    image_names = parse_image_names(images, N)
    @assert length(image_names) == length(x_avs) "$length(image_names) != $length(x_avs)"

    println("Building $(length(x_avs)) .mat files from 1 batch...")
    for i in eachindex(x_avs)
        mat_dict = Dict(
            "x" => x_avs[i],
            "y" => y_avs[i],
            "pass_sizes" => raw_piv_results[i][3],
            "overlap" => ol,
            "method" => "multin",
            "fn" => images,
            "u" => u_avs[i],
            "v" => v_avs[i],
            "npts" => npts[i],
            "uStd" => u_stds[i],
            "vStd" => v_stds[i]
        )
        # Crazy indexing to get the correct filename
        MAT.matwrite("$out_dir$(image_names[i]).mat", mat_dict)
    end


end

# MULTIBATCH
function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_dir::String, multi_batch::Int32) where {T}

    batches::Vector{String} = readdir(in_dir)
    mat_file_count = 0
    
    # Image pre-processing
    for batch in batches
        in_path = "$in_dir$batch" 
        images = get_raw_images(in_path)
        if length(images) ÷ N <= 1
            error("\n\nNumber of images in batch ($(length(images))) not divisible by N ($N).\n\n")
        end
        cropped_pairs = crop_and_pair_images(images, crop_factor)

        # Preallocate results from PIV
        raw_piv_results = Vector{Tuple{Tuple{Matrix{T}, Matrix{T}}, Tuple{Matrix{T}, 
        Matrix{T}}, Vector{Int32}} where {T}}()

        # Thread PIV processing
        Threads.@threads for pair in cropped_pairs
            push!(raw_piv_results, JuliaPIV.main(pair, Int32(final_win_size), Float32(ol)))
        end

        # Run stats
        ((x_avs, y_avs), (u_avs, v_avs), (u_stds, v_stds), npts) = statistics_of_piv(raw_piv_results, N)

        image_names = parse_image_names(images, N)
        @assert length(image_names) == length(x_avs) "$length(image_names) != $length(x_avs)"

        println("Building $(length(x_avs)) .mat files...")
        # Write to .mat file at argued out_dir
        for i in eachindex(x_avs)
            mat_dict = Dict(
            "x" => x_avs[i],
            "y" => y_avs[i],
            "pass_sizes" => raw_piv_results[i][3],
            "overlap" => ol,
            "method" => "multin",
            "fn" => images,
            "u" => u_avs[i],
            "v" => v_avs[i],
            "npts" => npts[i],
            "uStd" => u_stds[i],
            "vStd" => v_stds[i]
            )
            MAT.matwrite("$out_dir$(image_names[i]).mat", mat_dict)
            mat_file_count += 1
        end
    end
    println("Wrote $mat_file_count .mat files from $(length(batches)) batches.")
end

function parse_six_args()
        N = parse(Int32, ARGS[1])

        # Parse and split up crop factors
        crop_factors = split(ARGS[2], ",")
        crop_factors = [strip(crop_factor) for crop_factor in crop_factors]
        crop_factors = [parse.(Int32, crop_factor) for crop_factor in crop_factors]
        @assert length(crop_factors) == 4 "There should be 4 crop factors!"
        crop_factors = (crop_factors[1], crop_factors[2], crop_factors[3], crop_factors[4])

        final_win_size = parse(Int32, ARGS[3])
        ol = parse(Float32, ARGS[4])
        out_dir = "tests/pipeline_utility_testing/SVSout_23227179_1724441851/pivframes/"
        in_path = "tests/pipeline_utility_testing/testbatches/tmp.cp0atEmCuS.txt"
        return N, crop_factors, final_win_size, ol, out_dir, in_path
end

function parse_seven_args()
    N, crop_factors, final_win_size, ol, out_dir = parse_six_args()
    multi_batch = parse(Int32, ARGS[7])
    return N, crop_factors, final_win_size, ol, out_dir, ARGS[6], multi_batch
end

function julia_main()::Cint
    # Check on ARGS
    if length(ARGS) == 6
        N, crop_factors, final_win_size, ol, out_dir, in_path = parse_six_args()
        # Run PIV pipeline
        try 
            io_main(N, crop_factors, final_win_size, ol, out_dir, in_path)
        catch e
            error(e)
            return 1
        end
    elseif length(ARGS) == 7
        N, crop_factors, final_win_size, ol, out_dir, in_dir, multi_batch = parse_seven_args()
        # Run PIV pipeline
        try 
            io_main(N, crop_factors, final_win_size, ol, out_dir, in_dir, multi_batch)
        catch e
            error(e)
            return 1
        end
    else
        error("\nIncorrect number of arguments! Should be:\n-N\n-crop_factors\n-final_win_size\n-ol\n-out_dir\n-in_dir\n")
        return 1
    end
    return 0
end

julia_main()
end