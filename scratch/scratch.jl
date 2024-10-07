using CUDA
using StaticArrays

# function test_kernel!(U2, histo, histostd, tmp)
function test_kernel!(U2)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    j = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    
    tmp = @MArray Float32[NaN NaN NaN; 
                          NaN NaN NaN; 
                          NaN NaN NaN]

    # TODO: Find a way around slicing
    if i > 1 && i < size(U2, 1) && j > 1 && j < size(U2, 2)
        @inbounds begin
            tmp[1, 1] = U2[j-1, i-1]
            tmp[1, 2] = U2[j-1, i]
            tmp[1, 3] = U2[j-1, i+1]
            tmp[2, 1] = U2[j, i-1]
            tmp[2, 3] = U2[j, i+1]
            tmp[3, 1] = U2[j+1, i-1]
            tmp[3, 2] = U2[j+1, i]
            tmp[3, 3] = U2[j+1, i+1]
        end
    end
    return
end

# function bencher(U2, histo, histostd, tmp)
function bencher(U2)
    # 2D array size
    M = size(U2, 1)
    N = size(U2, 2)

    # Get config for kernel
    # kernel = @cuda launch=false test_kernel!(U2, histo, histostd, tmp)
    kernel = @cuda launch=false test_kernel!(U2)
    config = launch_configuration(kernel.fun)

    # Specificy number of threads and blocks for 2D array
    threads_x = min(M, config.threads)
    blocks_x = cld(M, threads_x)
    threads_y = min(N, config.threads)
    blocks_y = cld(N, threads_y)
    
    CUDA.@sync begin
        # kernel(U2, histo, histostd, tmp; 
        kernel(U2; 
               threads=(threads_x, threads_y), blocks=(blocks_x, blocks_y)
            )
    end

    return
end