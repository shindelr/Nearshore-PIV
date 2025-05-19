module PIVPipelineUtility

using FileIO
using Images
using Statistics
using MAT
include("./main.jl")

# Remove before compiling! -----------------------------------------------------
const ARGS = ["3", "1, 3072, 1, 2048", "16", "0.5", 
            "/home/server/pi/homes/shindelr/2023-test/piv-mat-out/custom-frame-verify/run4/",       # Output
            "/home/server/pi/homes/shindelr/2023-test/custom-frame-verify-jpgs/custom-frame-verify.txt", "1"]   # Input

"""
    get_raw_images(path::String)::Vector{String}

    Read in a list of image names from a given text file. Note that if 
    the text file is not in the same directory as the images, the path to 
    the images must be prepended to the image names.   

    Arguments:
        - `path::String`: Relative path to the text file containing image names.

    Returns:
        - `Vector{String}`: Vector of image names.
"""
function get_raw_images(path::String, N::Int32)::Vector{String}
    files::Vector{String} = readlines(path)
    @assert length(files) % N == 0 "Number of files in batch ($(length(files))) should be divisible by N ($(N))" 
    return files
end

"""
    crop_and_pair_images(images::Vector{String}, 
                        crop_factor::NTuple{4, Int32}
                        )::Vector{Tuple{Matrix{Float32}, Matrix{Float32}}}

    Create pairs of images from a given list of images. While pairing,
    the images are cropped and converted to Gray scale type.

    Arguments:
        - `images::Vector{String}`: Vector of image names.
        - `crop_factor::NTuple{4, Int32}`: Tuple of 4 integers representing the 
            cropping factor (left, right, top, bottom).
    Returns:
        - `Vector{Tuple{Matrix{Float32}, Matrix{Float32}}}`: A vector of
            tuples containing the processed images.
"""
function crop_and_pair_images(images::Vector{String}, crop_factor::NTuple{4, Int32}
                            )::Vector{Tuple{Matrix{Float32}, Matrix{Float32}}}
    # Preallocate image pair types
    image_pairs = Vector{Tuple{Matrix{Float32}, Matrix{Float32}}}()
    i = 1
    while i < length(images)
        # Load images and crop
        img1 = load(images[i])
        img2= load(images[i+1])
        img1 = img1[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        img2 = img2[crop_factor[3]:crop_factor[4], crop_factor[1]:crop_factor[2]]
        push!(image_pairs, (img1, img2))
        i += 2
    end
    @assert length(image_pairs) == length(images) ÷ 2 "Length of image pairs $(length(image_pairs)) should be half the length of images $(length(images))"
    return image_pairs
end

"""
    crop_and_group_images(images::Vector{String}, 
                        crop_factor::NTuple{4, Int32}, N::Int32
                        )::Vector{Vector{Matrix{Gray{N0f8}}}}

    Create non-pair groups of images from a given list of images. While 
    grouping, the images are cropped and converted to Gray scale type.
    This function is important for processing groups larger than 2.

    Arguments:
        - `images::Vector{String}`: Vector of image names.
        - `crop_factor::NTuple{4, Int32}`: Tuple of 4 integers representing the 
            cropping factor (left, right, top, bottom).
        - `N::Int32`: Number of images in each group.
    Returns:
        - `Vector{Vector{Matrix{GrayN0f8}}}}}`: A vector of
            vectors containing groups of the processed images.
"""
function crop_and_group_images(images::Vector{String}, 
                                crop_factor::NTuple{4, Int32},
                                N::Int32
                                )::Vector{Vector{Matrix{Gray{N0f8}}}}

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
            # Push images into their respective group
            # push!(group, Gray.(img)) 
            push!(group, img) 
            j += 1
        end
        # Push groups into the final return vector
        push!(image_groups, group)
        i += N
        count += 1
    end
    @assert length(image_groups) == length(images) ÷ N "Number of groups should be $(length(images)) ÷ $N"
    return image_groups
end

"""
    statistics_of_piv_groups(piv_results, N::Int32)

    Compute the statistics of PIV results for variable groups of
    images. After partitioning the results of the PIV algorithm into
    subgroups, each group's data is averaged and the standard deviation
    calculated.

    Arguments:
        - `piv_results`: An array containing each of the piv results 
            coming back from the main PIV algorithm.
    Returns:
        - `Tuple{Tuple{}, Tuple{}, Matrix{Int32}}`: A tuple of 
            tuples containing the averages, standard deviations, and a 
            sum of points that were NaN prior to averaging for each 
            pair of images.
"""
function statistics_of_piv_groups(piv_results, N::Int32)
    println("Calculating statistics...")
    # Preallocations
    u_avs = Vector{Matrix{Float32}}(); u_stds = Vector{Matrix{Float32}}()
    v_avs = Vector{Matrix{Float32}}(); v_stds = Vector{Matrix{Float32}}()
    x_choices = Vector{Matrix{Float32}}(); y_choices = Vector{Matrix{Float32}}() 
    npts = Vector{Matrix{Int32}}()

    # Unpack results
    us, vs, = Vector{Matrix{Float32}}(), Vector{Matrix{Float32}}()
    xs, ys = Vector{Matrix{Float32}}(), Vector{Matrix{Float32}}()
    for result in piv_results
        push!(xs, result[1][1])
        push!(ys, result[1][2])
        push!(us, result[2][1])
        push!(vs, result[2][2])
    end

    
    group = 1
    i = 1
    offset = length(piv_results) ÷ N
    while group <= N
        u_group = us[i:(i + offset - 1)]
        v_group = vs[i: (i + offset - 1)]
        # Really nifty bitmask summation of NaN locations in each matrix
        nans = reduce(+, [isnan.(u) for u in u_group])
        push!(u_avs, nan_mean(u_group))
        push!(u_stds, nan_std(u_group))
        push!(v_avs, nan_mean(v_group))
        push!(v_stds, nan_std(v_group))
        push!(npts, nans)
        push!(x_choices, xs[i])
        push!(y_choices, ys[i])
        i += offset
        group += 1
    end

    return ((x_choices, y_choices), (u_avs, v_avs), (u_stds, v_stds), npts)
end

"""
    nan_mean(arr::Vector{Matrix{Float32}})::Matrix{Float32}

    Compute the mean of a vector of matrices, ignoring NaN values. 

    Arguments:
        - `arr::Vector{Matrix{Float32}}`: Vector of matrices containing 
            the data to be averaged. May contain NaN values.
    Returns:
        - `Matrix{Float32}`: A single matrix containing the mean of all the 
            input matrices disregarding NaN values.
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

    Compute the standard deviation of a vector of matrices, ignoring 
    NaN values. NOTE: If there is only one value in the vector, the
    standard deviation is undefined and set to 0.0. If all values in
    the matrix return NaN, the function will halt and throw an error.

    Arguments:
        - `arr::Vector{Matrix{Float32}}`: Vector of matrices containing 
            the data to be averaged. May contain NaN values.
    Returns:
        - `Matrix{Float32}`: A single matrix containing the std of all the 
            input matrices disregarding NaN values.
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
            elseif length(temp_stds) == 1
                std_val = 0.0  # If only one value, std = 0 bc no variance for one value
            end
            std_matrix[i, j] = std_val
        end
    end
    @assert all(isnan.(std_matrix)) == false "All values in std_matrix should not be NaN"
    return std_matrix
end

"""
    parse_image_names(images::Vector{String}, N::Int32)::Vector{Vector{String}}

    Parse a vector of image names into groups of N image names. 

    Arguments:
        - `images::Vector{String}`: Vector of image names to be 
            parsed into groups
        - `N::Int32`: Number of images expected to be in each group.
    Returns:
        - `Vector{Vector{String}}`: A vector of groups of image names.
"""
function parse_image_names(images::Vector{String}, N::Int32)::Vector{Vector{String}}
    # Take just the name of the image itself and not full filepath
    image_names = [basename(image) for image in images]
    image_groups = Vector{Vector{String}}()
    for group in Iterators.partition(image_names, N)
        push!(image_groups, group)
    end
    return image_groups
end

"""

    Run PIV on a batch of images where N=2.

"""
function paired_piv(N::T, final_win_size::T, ol::Float32, out_dir::String, 
                    images::Vector{String}, crop_factor::Tuple{T,T,T,T}) where {T}

    # Preallocate results from PIV: [(x, y), (u, v), pass_sizes]
    raw_piv_results = Vector{Tuple{
                            Tuple{Matrix{T}, Matrix{T}}, 
                            Tuple{Matrix{T}, Matrix{T}}, 
                            Vector{Int32}
                            } where {T}}()

    cropped_pairs = crop_and_pair_images(images, crop_factor)
    @assert length(cropped_pairs) == length(images) ÷ 2 "Length of cropped pairs should be half the length of images"
    image_groups_names = parse_image_names(images, N)

    println("Running PIV...")
    for pair in cropped_pairs
        # Run PIV proper!
        push!(raw_piv_results, main(pair, Int32(final_win_size), Float32(ol)))
    end

    println("Building $(length(raw_piv_results)) .mat files from 1 batch...")
    pass_sizes = [raw_piv_results[1][3] raw_piv_results[1][3]]  # Just a formatting thing to match OG Matlab
    for (i, result) in enumerate(raw_piv_results)
        group_name = image_groups_names[i]
        x = result[1][1]
        y = result[1][2]
        u = result[2][1]
        v = result[2][2]
        npts = isnan.(u)
        mat_dict = Dict(
            "x" => x,
            "y" => y,
            "pass_sizes" => pass_sizes,
            "overlap" => ol,
            "method" => "multin",
            "fn" => group_name,
            "u" => u,
            "v" => v,
            "npts" => npts,
        )
        MAT.matwrite("$out_dir/$(group_name[1]).mat", mat_dict)
    end
end

"""
    Run PIV when N > 2.
"""
function grouped_piv(N::T, final_win_size::T, ol::Float32, out_dir::String, 
                    images::Vector{String}, crop_factor::Tuple{T,T,T,T}) where {T}

    # Preallocate results from PIV: [(x, y), (u, v), pass_sizes]
    raw_piv_results = Vector{Tuple{
                             Tuple{Matrix{T}, Matrix{T}}, 
                             Tuple{Matrix{T}, Matrix{T}}, 
                             Vector{Int32}
                             } where {T}}()

    image_groups = crop_and_group_images(images, crop_factor, N)
    for group in image_groups
        for i in (1:length(group)-1)
            push!(raw_piv_results, main((group[i], group[i+1]), Int32(final_win_size), Float32(ol)))
        end
    end

    # Explicitly setting subgroup size for clarity
    subgroup_size = Int32(length(image_groups))
    
    # PIV stats
    ((xs, ys),
    (u_avs, v_avs),
    (u_stds, v_stds), 
    npts) = statistics_of_piv_groups(raw_piv_results, subgroup_size)

    # Format pass_sizes and group image names for .mat file
    image_groups_names = parse_image_names(images, N)
    @assert length(image_groups_names) == length(u_avs) "$(length(image_groups_names)) != $(length(u_avs))"
    
    println("Building $(length(u_avs)) .mat files from 1 batch...")
    pass_sizes = [raw_piv_results[1][3] raw_piv_results[1][3]]
    for i in eachindex(u_avs)
        mat_dict = Dict(
            "x" => xs[i],
            "y" => ys[i],
            "pass_sizes" => pass_sizes,
            "overlap" => ol,
            "method" => "multin",
            "fn" => image_groups_names[i],  
            "u" => u_avs[i],
            "v" => v_avs[i],
            "npts" => npts[i],
            "uStd" => u_stds[i],
            "vStd" => v_stds[i]
        )
        MAT.matwrite("$out_dir/$(image_groups_names[i][1]).mat", mat_dict)
    end

end

"""
    io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String) where {T}

    Run a single batch of images through the Julia PIV algorithm.

    Arguments:
        - `N::Int32`: Number of images in each subgroup to run PIV on.
            Corresponds to the LiDAR scan rate.
        - `crop_factor::Int32`: Tuple of 4 integers representing the 
            cropping factor (left, right, top, bottom).
        - `final_win_size::Int32`: Final window size for PIV.
        - `ol::Float32`: Overlap percentage for PIV.
        - `out_dir::String`: Directory to write .mat files to.
        - `in_path::String`: Path to the directory containing the images.

    Returns:
        None
"""
function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float32, out_dir::String, in_path::String) where {T}

    # Image pre-processing
    images = get_raw_images(in_path, N)
    if length(images) ÷ N <= 1
        error("\n\nNumber of images in directory ($(length(images))) not divisible by N ($N).\n\n")
    end

    if N == 2
        paired_piv(N, final_win_size, ol, out_dir, images, crop_factor)
        return 
    elseif N > 2
        grouped_piv(N, final_win_size, ol, out_dir, images, crop_factor)
        return
    else
        error("N should be greater than 1")
    end
end

"""
    parse_seven_args()

    Parse arguments from command line.

    Returns:
        - `Tuple{Int32, NTuple{4, Int32}, Int32, 
                Float32, String, String, Int32}`: ARGS to run JuliaPIV.
"""
function parse_seven_args()
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
        verbose = parse(Int32, ARGS[7])
        return N, crop_factors, final_win_size, ol, out_dir, in_path, verbose
end


"""
    julia_main()::Cint

    Main entry point to run the PIV pipeline utility. Two optional flags
    can be passed in from the command line. `verbose` will suppress all
    stdout if set to 0. `multi_batch` will run multiple batches of images
    through the PIV algorithm if set to 1. Please note, that when running
    multiple batches, the `in_path` argument should be the directory
    containing the batches, not the batch itself.

    .mat files will be written
    to the argued `out_dir` directory path. These .mat files will contain
    a variety of information detailed here:
        x: [255×299 double]
        y: [255×299 double]
        pass_sizes: [3×2 double]
        overlap: 0.5
            method: 'multin'
            fn: {list of jpg files}
                u: [255×299 double]
                v: [255×299 double]
            npts: [255×299 double]  # number of data points that weren't NaN 
                                    # prior to time-average
            uStd: [255×299 double]  # standard deviation of the N results
            vStd: [255×299 double]  # ditto

    Returns:
        - `Cint`: 0 if successful, 1 if unsuccessful.

"""
function julia_main()::Cint
    # Check on ARGS
    if length(ARGS) == 7
        N, crop_factors, final_win_size, ol, out_dir, in_path, verbose = parse_seven_args()
        if verbose == 0
            og_stdout = stdout
            redirect_stdout(devnull)
        end
        # Run PIV pipeline
        try 
            io_main(N, crop_factors, final_win_size, ol, out_dir, in_path)
        catch e
            error(e)
            return 1
        end
    else
        println(ARGS)
        error("\nIncorrect number of arguments! Should be:\n-N\n-crop_factors\n-final_win_size\n-ol\n-out_dir\n-in_dir\n")
        return 1
    end
    if verbose == 0
        redirect_stdout(og_stdout)
    end
    return 0
end

end