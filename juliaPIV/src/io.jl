module PIVPipelineUtility
export julia_main 

using Base.Threads
using FileIO
using Images
using Plots           # Can be removed
using Statistics
using MAT
using DelimitedFiles  # Can be removed
include("main.jl")
using .JuliaPIV

"""
    get_raw_images(PATH::String)::Vector{String}

    Briefly describe what the function does.

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function get_raw_images(PATH::String)::Vector{String}
    files::Vector{String} = readlines(PATH)
    prefix_dir = "tests/pipeline_utility_testing/"
    # Get raw images and prepend the test directory
    return ["$prefix_dir$file" for file in files]
end

"""
    crop_and_pair_images(images::Vector{String}, 
                        crop_factor::NTuple{4, Int32}
                        )::Vector{Tuple{Matrix{Float32}, Matrix{Float32}}}

    Briefly describe what the function does.

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function crop_and_pair_images(images::Vector{String}, crop_factor::NTuple{4, Int32})::Vector{Tuple{Matrix{Float32}, Matrix{Float32}}}
    image_pairs = Dict{Int32, Tuple{String, String, Matrix{Gray{N0f8}}, Matrix{Gray{N0f8}}}}()
    i = 1
    count = 1
    while i < length(images)
        img1 = load(images[i])
        img2= load(images[i+1])
        img1 = img1[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        img2 = img2[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        image_pairs[count] = (images[i], images[i+1], Gray.(img1), Gray.(img2))
        i += 2
        count += 1
    end

    @assert length(image_pairs) == length(images) ÷ 2 "Length of image pairs should be half the length of images"
    return [(image_pairs[i][3], image_pairs[i][4]) for i in eachindex(image_pairs)]
end

"""
    crop_and_group_images(images::Vector{String}, 
                        crop_factor::NTuple{4, Int32}, N::Int32
                        )::Vector{Vector{Matrix{Gray{N0f8}}}}

    Briefly describe what the function does.

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function crop_and_group_images(images::Vector{String}, crop_factor::NTuple{4, Int32}, N::Int32)
    # image_groups = Dict{Int32, Vector{Matrix}}()
    image_groups = Vector{Vector{Matrix{Gray{N0f8}}}}()
    i = 1
    count = 1
    while i < length(images)
        j = 1
        group = Vector{Matrix{Gray{N0f8}}}()
        # group --> Vector(matrix_j)
        while j <= N
            img_name = images[i + j - 1]
            img = load(img_name)
            img = img[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
            push!(group, Gray.(img)) 
            j += 1
        end
        # image_groups[count] = group
        push!(image_groups, group)
        i += N
        count += 1
    end
    @assert length(image_groups) == length(images) ÷ N "Number of groups should be $(length(images)) ÷ $N"
    return image_groups
end

"""
    statistics_of_piv_pairs(piv_results)

    Briefly describe what the function does.

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function statistics_of_piv_pairs(piv_results)
    println("Calculating statistics...")
    # Preallocations
    u_avs = Vector{Matrix{Float32}}()
    u_stds = Vector{Matrix{Float32}}()
    v_avs = Vector{Matrix{Float32}}()
    v_stds = Vector{Matrix{Float32}}()
    x_avs = Vector{Matrix{Float32}}()
    y_avs = Vector{Matrix{Float32}}()
    npts = Vector{Matrix{Int32}}()

    # Unpack results
    xs, ys, us, vs = [], [], [], []
    for result in piv_results
        push!(xs, result[1][1])
        push!(ys, result[1][2])
        push!(us, result[2][1])
        push!(vs, result[2][2])
    end

    # Calculate statistics
    i = 1
    while i <= length(us)
        # Make groups of 2 results for averages, stds, etc.
        un_group = us[i:min(i+1, length(us))]
        vn_group = vs[i:min(i+1, length(us))]
        xn_group = xs[i:min(i+1, length(us))]
        yn_group = ys[i:min(i+1, length(us))]
        
        # Create binary masks for each matrix in u-group
        nan_binary_masks = Vector{Matrix{Bool}}([isnan.(u) for u in un_group])
        push!(npts, sum(nan_binary_masks))

        # Compute averages and stds, no need to worry about NaNs for pairs
        push!(u_avs, mean(un_group))
        push!(u_stds, std(un_group))
        push!(v_avs, mean(vn_group))
        push!(v_stds, std(vn_group))
        push!(x_avs, mean(xn_group))
        push!(y_avs, mean(yn_group))

        i += 1
    end
    
    return ((x_avs, y_avs), (u_avs, v_avs), (u_stds, v_stds), npts)
end

"""
    statistics_of_piv_groups(piv_results, N::Int32)

    Briefly describe what the function does.

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function statistics_of_piv_groups(piv_results, N::Int32)
    println("Calculating statistics...")
    # Preallocations
    u_avs = Vector{Matrix{Float32}}()
    u_stds = Vector{Matrix{Float32}}()
    v_avs = Vector{Matrix{Float32}}()
    v_stds = Vector{Matrix{Float32}}()
    x_avs = Vector{Matrix{Float32}}()
    y_avs = Vector{Matrix{Float32}}()
    npts = Vector{Matrix{Int32}}()

    # Unpack results
    xs, ys, us, vs = [], [], [], []
    for result in piv_results
        push!(xs, result[1][1])
        push!(ys, result[1][2])
        push!(us, result[2][1])
        push!(vs, result[2][2])
    end

    # Calculate statistics, partition groups by size N
    for (un_group::Vector{Matrix{Float32}}, vn_group::Vector{Matrix{Float32}}, 
        xn_group::Vector{Matrix{Float32}}, yn_group::Vector{Matrix{Float32}}) in zip(
        Iterators.partition(us, N), 
        Iterators.partition(vs, N), 
        Iterators.partition(xs, N), 
        Iterators.partition(ys, N))

        # Create binary masks for each matrix in u-group
        nan_binary_masks = Vector{Matrix{Bool}}([isnan.(u) for u in un_group])
        push!(npts, sum(nan_binary_masks))

        # Compute averages and stds
        push!(u_avs, nan_mean(un_group))
        push!(u_stds, nan_std(un_group))
        push!(v_avs, nan_mean(vn_group))
        push!(v_stds, nan_std(vn_group))
        push!(x_avs, nan_mean(xn_group))
        push!(y_avs, nan_mean(yn_group))
    end
    
    return ((x_avs, y_avs), (u_avs, v_avs), (u_stds, v_stds), npts)
end

"""
    nan_mean(arr::Vector{Matrix{Float32}})::Matrix{Float32}

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function nan_mean(arr::Vector{Matrix{Float32}})::Matrix{Float32}
    # Preallocate
    mean_matrix = Matrix{Float32}(undef, size(arr[1]))
    for i in 1:size(arr[1], 1)
        for j in 1:size(arr[1], 2)
            # If all values are NaN, mean_val = NaN
            mean_val = NaN
            temp_means = Vector{Float32}()
            for m in arr
                if !isnan(m[i, j])  # Only want the mean of non-nans
                    push!(temp_means, m[i, j])
                end
            end
            if length(temp_means) > 0
                mean_val = mean(temp_means)
            end
            mean_matrix[i, j] = mean_val
        end
    end
    @assert all(isnan.(mean_matrix)) == false "All values in mean_matrix should not be NaN"
    return mean_matrix
end

"""
    nan_std(arr::Vector{Matrix{Float32}})::Matrix{Float32}

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function nan_std(arr::Vector{Matrix{Float32}})::Matrix{Float32}
    std_matrix = Matrix{Float32}(undef, size(arr[1]))
    for i in 1:size(arr[1], 1)
        for j in 1:size(arr[1], 2)
            std_val = NaN
            temp_stds = Vector{Float32}()
            for m in arr
                if !isnan(m[i, j])
                    push!(temp_stds, m[i, j])
                end
            end
            if length(temp_stds) > 1
                std_val = std(temp_stds)
            # elseif length(temp_stds) == 1
            #     std_val = 0.0  # If only one value, std = 0 bc no variance for one value
            end
            std_matrix[i, j] = std_val
        end
    end
    @assert all(isnan.(std_matrix)) == false "All values in std_matrix should not be NaN"
    return std_matrix
end

"""
    parse_image_names(images::Vector{String}, N::Int32)::Vector{String}

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
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
"""
    io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String) where {T}

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String) where {T}

    # Preallocate results from PIV
    raw_piv_results = Vector{Tuple{
                             Tuple{Matrix{T}, Matrix{T}}, 
                             Tuple{Matrix{T}, Matrix{T}}, 
                             Vector{Int32}
                             } where {T}}()

    # Image pre-processing
    images = get_raw_images(in_path)
    if length(images) ÷ N <= 1
        error("\n\nNumber of images in directory ($(length(images))) not divisible by N ($N).\n\n")
    end

    if N == 2
        cropped_pairs = crop_and_pair_images(images, crop_factor)
        @assert length(cropped_pairs) == length(images) ÷ 2 "Length of cropped pairs should be half the length of images"
        Threads.@threads for pair in cropped_pairs
            push!(raw_piv_results, JuliaPIV.main(pair, Int32(final_win_size), Float32(ol)))
        end
        # PIV stats
        ((x_avs, y_avs),
        (u_avs, v_avs), 
        (u_stds, v_stds), 
         npts) = statistics_of_piv_pairs(raw_piv_results)

    elseif N > 2
        image_groups = crop_and_group_images(images, crop_factor, N)
        Threads.@threads for group in image_groups
            for i in (1:length(group)-1)
                push!(raw_piv_results, JuliaPIV.main((group[i], group[i+1]), Int32(final_win_size), Float32(ol)))
            end
        end

        # Explicitly getting subgroup size for clarity
        subgroup_size = Int32(length(image_groups))
        # PIV stats
        ((x_avs, y_avs), 
        (u_avs, v_avs), 
        (u_stds, v_stds), 
        npts) = statistics_of_piv_groups(raw_piv_results, subgroup_size)

    else
        error("N should be greater than 1")
    end

    image_names = parse_image_names(images, N)
    @assert length(image_names) == length(x_avs) "$(length(image_names)) != $(length(x_avs))"
    pass_sizes = [raw_piv_results[1][3] raw_piv_results[1][3]]

    println("Building $(length(x_avs)) .mat files from 1 batch...")
    for i in eachindex(x_avs)
        mat_dict = Dict(
            "x" => x_avs[i],
            "y" => y_avs[i],
            "pass_sizes" => pass_sizes,
            "overlap" => ol,
            "method" => "multin",
            "fn" => image_names,  
            "u" => u_avs[i],
            "v" => v_avs[i],
            "npts" => npts[i],
            "uStd" => u_stds[i],
            "vStd" => v_stds[i]
        )
        MAT.matwrite("$out_dir$(image_names[i]).mat", mat_dict)
    end
end

# MULTIBATCH
"""
    io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_dir::String, multi_batch::Int32) where {T}

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_dir::String, multi_batch::Int32) where {T}

    batches::Vector{String} = readdir(in_dir)
    mat_file_count = 0
    
    # Image pre-processing
    for batch in batches
        # Preallocate results from PIV
        raw_piv_results = Vector{Tuple{
                                Tuple{Matrix{T}, Matrix{T}}, 
                                Tuple{Matrix{T}, Matrix{T}}, 
                                Vector{Int32}
                                } where {T}}()
                  
        # Get batch of raw images
        in_path = "$in_dir$batch" 
        images = get_raw_images(in_path)
        if length(images) ÷ N <= 1
            error("\n\nNumber of images in batch ($(length(images))) not divisible by N ($N).\n\n")
        end

        if N == 2
            cropped_pairs = crop_and_pair_images(images, crop_factor)
            @assert length(cropped_pairs) == length(images) ÷ 2 "Length of cropped pairs should be half the length of images"
            Threads.@threads for pair in cropped_pairs
                push!(raw_piv_results, JuliaPIV.main(pair, Int32(final_win_size), Float32(ol)))
            end
            # PIV stats
            ((x_avs, y_avs),
            (u_avs, v_avs), 
            (u_stds, v_stds), 
             npts) = statistics_of_piv_pairs(raw_piv_results)
    
        elseif N > 2
            image_groups = crop_and_group_images(images, crop_factor, N)
            Threads.@threads for group in image_groups
                for i in (1:length(group)-1)
                    push!(raw_piv_results, JuliaPIV.main((group[i], group[i+1]), Int32(final_win_size), Float32(ol)))
                end
            end
            # Explicitly getting subgroup size for clarity
            subgroup_size = Int32(length(image_groups))
            # PIV stats
            ((x_avs, y_avs), 
            (u_avs, v_avs), 
            (u_stds, v_stds), 
            npts) = statistics_of_piv_groups(raw_piv_results, subgroup_size)

        else
            error("N should be greater than 1")
        end

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
            "fn" => image_names,
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

"""
    parse_six_args()

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
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
        out_dir = ARGS[5]
        in_path = ARGS[6]
        return N, crop_factors, final_win_size, ol, out_dir, in_path
end

"""
    parse_seven_args()

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
function parse_seven_args()
    N, crop_factors, final_win_size, ol, out_dir = parse_six_args()
    multi_batch = parse(Int32, ARGS[7])
    return N, crop_factors, final_win_size, ol, out_dir, ARGS[6], multi_batch
end

"""
    julia_main()::Cint

    Arguments:
        - `arg1::Type`: Description of the first argument.
        - `arg2::Type`: Description of the second argument.
        - `args...`: Description of additional arguments.

    Returns:
        - `ReturnType`: Description of the return value.
"""
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