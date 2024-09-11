using Base.Threads
using FileIO
using Images
using Plots
using Statistics
using MAT
include("main.jl")
using .JuliaPIV

function get_raw_images(PATH::String)::Vector{String}
    files::Vector{String} = readlines(PATH)
    prefix_dir = "tests/pipeline_utility_testing/"
    # Get raw images and prepend the test directory
    return ["$prefix_dir$file" for file in files]
end

function crop_and_pair_images(images::Vector{String}, crop_factor::Tuple{Int,Int,Int,Int})
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

function statistics_of_piv(piv_results, n::Int)
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
    t_avs = Vector{Matrix{Float32}}()

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
        un_group = us[i:min(i+n-1, length(us))]
        vn_group = vs[i:min(i+n-1, length(us))]
        xn_group = xs[i:min(i+n-1, length(us))]
        yn_group = ys[i:min(i+n-1, length(us))]
        push!(u_avs, mean(un_group))
        push!(u_stds, std(un_group))
        push!(v_avs, mean(vn_group))
        push!(v_stds, std(vn_group))
        push!(x_avs, mean(xn_group))
        push!(t_avs, mean(yn_group))
        i += n
    end

    return ((x_avs, t_avs), (u_avs, v_avs), (u_stds, v_stds))
end

function io_main(N::T, crop_factor::Tuple{T,T,T,T}, final_win_size::T,
    ol::Float64, out_dir::String, in_path::String) where {T}

    # Image pre-processing
    images = get_raw_images(in_path)
    cropped_pairs = crop_and_pair_images(images, crop_factor)

    # PIV, could thread here?
    raw_piv_results = Vector{Tuple{Tuple{Matrix{T}, Matrix{T}}, Tuple{Matrix{T}, Matrix{T}}, Vector{Int32}} where {T}}()

    Threads.@threads for pair in cropped_pairs
        push!(raw_piv_results, JuliaPIV.main(pair, Int32(final_win_size), Float32(ol)))
    end

    # (npts, u_std, v_std, u_av, v_av) = statistics_of_piv(raw_piv_results)
    ((x_avs, y_avs), (u_avs, v_avs), (u_stds, v_stds)) = statistics_of_piv(raw_piv_results, N)

    i = 1
    image_names = Vector{String}()
    while i < length(images)
        push!(image_names, images[i][69:end-4])
        i += N
    end
    @assert length(image_names) == length(x_avs)

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
            "npts" => "TODO",
            "uStd" => u_stds[i],
            "vStd" => v_stds[i]
        )
        # Crazy indexing to get the correct filename
        MAT.matwrite("$out_dir$(image_names[i]).mat", mat_dict)
    end


end


# stats_piv_final_res = {Tuple(
#                   Tuple(X,Y),
#                   Tuple(U,V),
#                   Vec(pass_sizes),
#                   Vec(jpg files),
#                   Matrix(nan-points),
#                   Matrix(uStd),
#                   Matrix(vStd),   
#                ),
#            Tuple(...),
#            ...
#            }

# Testing:
in_args = (
    N=2,
    crop_factor=(24, 2424, 1, 2048),
    final_win_size=16,
    ol=0.5,
    out_dir="tests/pipeline_utility_testing/SVSout_23227179_1724441851/pivframes/",
    in_path="tests/pipeline_utility_testing/testbatches/tmp.cp0atEmCuS.txt"
    # in_dir = readdir("tests/pipeline_utility_testing/testbatches")
)

timed() = @time io_main(in_args[1], in_args[2], in_args[3], in_args[4], 
                            in_args[5], in_args[6])
timed()