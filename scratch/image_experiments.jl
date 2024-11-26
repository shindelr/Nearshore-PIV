"""
Working on a pass kernel.

NOTES:
- Does not account for overlap at this point.
"""

using CUDA
using BenchmarkTools
using Images
using FileIO


const PIVWIN::Int32 = 16
const CROPS::NTuple{4, Int32} = (24, 2424, 1, 2048)

function main()
    # Load and crop image, allocate on GPU
    h_im_a = load("../juliaPIV/data/im1.jpg")
    h_im_a::Matrix{Float32} = convert(Matrix{Float32}, h_im_a)
    crop!(h_im_a, CROPS)
    d_im_a::CuArray{Float32} = cu(h_im_a)
    d_im_a_out::CuArray{Float32} = CUDA.zeros(size(d_im_a))
    m::Int32, n::Int32 = size(h_im_a)

    display(h_im_a)

    # Get sizes of everything for the kernel
    pass_sizes = get_passes()
    threads = (pass_sizes[1], pass_sizes[1])  # will need to iterate this

    # Janky tricks for the 64x64 blocks, which are too large for the individual
    # SMs on the GPU (max 32x32). First pass will need to use a kernel which
    # iterates over a 2x2 tile within the block per thread. So, the kernel still
    # needs the same size grid dimensions.
    if threads[1] == 64
        threads = Int32.(threads ./ 2)
        dim_grid = (Int32(ceil(m/threads[1])), Int32(ceil(n/threads[2])))
        dim_grid = Int32.(dim_grid ./ 2)
    end
    @cuda threads=threads blocks=dim_grid shmem=sizeof(d_im_a) _block_means(d_im_a_out, d_im_a, m, n)

    test_mean_block_kernel!(h_im_a)

    display(h_im_a)
    display(d_im_a_out)

    return isapprox(h_im_a, Array(d_im_a_out))
end

function _block_means(im_out::CuDeviceMatrix{Float32}, 
                    d_im_a::CuDeviceMatrix{Float32}, 
                    m::Int32, n::Int32)
    # Allocated block-wise shared memory for the sum
    block_sum = CUDA.@cuStaticSharedMem(Float32, 1)

    # Calculate global row, col in the image
    row = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    col = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    # Since each thread needs to process a 2x2 tile within the 32x32 block, we
    # accumulate a sum local to the tile before adding it to the block sum
    thread_local_sum = 0f0
    for i in 1:2, j in 1:2
        thread_row = Int32(row * 2 + j)
        thread_col = Int32(col * 2 + i)
        if thread_row <= m && thread_col <= n
           @inbounds thread_local_sum += d_im_a[thread_row, thread_col]
        end
    end

    # As threads finish with their local tile summations, we can atomically add
    # them to the block sum, which is a element vector, then compute the block
    # mean
    CUDA.@atomic block_sum[1] += thread_local_sum
    CUDA.sync_threads()  # Need to wait on all threads before computing mean 

    # Finally, scale each element in each thread's tile by the block mean
    for i in 1:2, j in 1:2
        thread_row = row * 2 + j
        thread_col = col * 2 + i
        if thread_row <= m && thread_col <= n
            @inbounds im_out[thread_row, thread_col] = d_im_a[thread_row, thread_col] - (block_sum[1] / (64^2))
        end
    end
    return nothing
end

function test_mean_block_kernel!(image) 
    for (cj, jj) in enumerate(1:64:1985), (ci, ii) in enumerate(1:64:3009)
        C = image[floor(Int32, jj):floor(Int32, jj + 64 - 1),
                floor(Int32, ii):floor(Int32, ii + 64 - 1)]
        
        C = C .- mean(C)
        image[floor(Int32, jj):floor(Int32, jj + 64 - 1), floor(Int32, ii):floor(Int32, ii + 64 - 1)] = C
    end
end

function crop!(image::Matrix{Float32}, crops::NTuple{4, Int32})
    image = image[crops[3]:crops[4], crops[1]:crops[2]]
    return nothing
end

function get_passes()
    log2pivwin::Int32 = log2(PIVWIN)
    pass_sizes::Vector{Int32} = 2 .^ collect(Int32, 6:-1:log2pivwin)
    push!(pass_sizes, PIVWIN)
    return pass_sizes
end

