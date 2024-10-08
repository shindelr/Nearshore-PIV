using CUDA
using StaticArrays

# function test_kernel!(U2, histo, histostd)
function test_kernel!(U2)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    j = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    
    subarray = @MVector fill(ComplexF32(NaN + NaN * im), 9)

    if i > 1 && i < size(U2, 1) && j > 1 && j < size(U2, 2)

        # Get subarray in a really ugly way
        @inbounds begin
            subarray[1] = U2[j-1, i-1]
            subarray[2] = U2[j-1, i]
            subarray[3] = U2[j-1, i+1]
            subarray[4] = U2[j, i-1]
            subarray[6] = U2[j, i+1]
            subarray[7] = U2[j+1, i-1]
            subarray[8] = U2[j+1, i]
            subarray[9] = U2[j+1, i+1]
        end

        # TODO: Compute im_median_magnitude and im_std

        # 1) Filter out NaN values
        not_nan_values = @MVector fill(ComplexF32(NaN + NaN * im), 9)
        not_nan = 0
        for index in 1:length(subarray)
            if !isnan(real(subarray[index])) && !isnan(imag(subarray[index])) 
                not_nan += 1
                not_nan_values[not_nan] = subarray[index]
            end
        end
        if not_nan == 0
            # Make histo[j, i] = NaN
            return
        end

        # 2) Sort based on abs2(x), then angle(x). Insertion sort algo

        # TODO: Implement insertion sort, having a bounds issue!
        for k in 2:not_nan
            l = k
            while l > 1 && abs2(not_nan_values[l]) > abs2(not_nan_values[l-1]) ||
                         (abs2(not_nan_values[l]) == abs2(not_nan_values[l-1]) &&
                          angle(not_nan_values[l]) > angle(not_nan_values[l-1]))
                
                temp = not_nan_values[l]
                not_nan_values[l] = not_nan_values[l - 1]
                not_nan_values[l - 1] = temp
                l -= 1
            end
        end

    end
    return
end

# function bencher(U2, histo, histostd)
function bencher()

    # U2::CuArray{ComplexF32} = cu([
    #        3 + 4im 1 + 2im 7 + 8im 0 + 0im 5 + 6im;
    #        9 + 10im 11 - 12im 15 - 16im 13 + 14im 7 - 8im;
    #        21 + 22im 17 + 18im 19 - 20im 25 + 26im 15 - 16im;
    #        27 - 28im 23 - 24im 29 + 30im 31 - 32im 25 + 26im
    #    ])
    U2 = @SMatrix ComplexF32[
           3 + 4im 1 + 2im 7 + 8im 0 + 0im 5 + 6im;
           9 + 10im 11 - 12im 15 - 16im 13 + 14im 7 - 8im;
           21 + 22im 17 + 18im 19 - 20im 25 + 26im 15 - 16im;
           27 - 28im 23 - 24im 29 + 30im 31 - 32im 25 + 26im
       ]

    # 2D array size
    M = size(U2, 1)
    N = size(U2, 2)


    # Get config for kernel
    # kernel = @cuda launch=false test_kernel!(U2, histo, histostd)
    kernel = @cuda launch=false test_kernel!(U2)
    config = launch_configuration(kernel.fun)

    # Specificy nukber of threads and blocks for 2D array
    threads_x = min(M, config.threads)
    blocks_x = cld(M, threads_x)
    threads_y = min(N, config.threads)
    blocks_y = cld(N, threads_y)
    
    CUDA.@sync begin
        # kernel(U2, histo, histostd; 
        kernel(U2; 
               threads=(threads_x, threads_y), blocks=(blocks_x, blocks_y)
            )
    end

    return
end