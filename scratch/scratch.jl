using CUDA
using StaticArrays

# function test_kernel!(U2, histo, histostd)
function test_kernel!(U2)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    j = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y
    
    subarray = @MVector fill(ComplexF32(NaN + NaN * im), 9)

    if i > 1 && i < size(U2, 1) && j > 1 && j < size(U2, 2)

        # 1) Get subarray in a really ugly way
        @inbounds begin
            subarray[1] = U2[j-1, i-1]
            subarray[2] = U2[j-1, i]
            subarray[3] = U2[j-1, i+1]
            subarray[4] = U2[j, i-1]
            subarray[5] = NaN + NaN * im
            subarray[6] = U2[j, i+1]
            subarray[7] = U2[j+1, i-1]
            subarray[8] = U2[j+1, i]
            subarray[9] = U2[j+1, i+1]
        end

        # 2) Filter out NaN values
        valid_vals = @MVector fill(ComplexF32(NaN + NaN * im), 9)
        not_nan_count = 0
        for index in 1:length(subarray)
            if !isnan(real(subarray[index])) && !isnan(imag(subarray[index])) 
                not_nan_count += 1
                valid_vals[not_nan_count] = subarray[index]
            end
        end
        if not_nan_count == 0
            # TODO: Make histo[j, i] = NaN
            return
        end

        # 3) Sort based on abs2(x), then angle(x). Insertion sort algo
        for k in 2:not_nan_count
            l = k
            while l > 1 && (abs2(valid_vals[l]) > abs2(valid_vals[l-1]) ||
                            (abs2(valid_vals[l]) == abs2(valid_vals[l-1]) &&
                             angle(valid_vals[l]) > angle(valid_vals[l-1])))
                
                temp = valid_vals[l]
                valid_vals[l] = valid_vals[l - 1]
                valid_vals[l - 1] = temp
                l -= 1
            end
        end

        if threadIdx().x == 2 && threadIdx().y == 2 && blockIdx().x == 1 && blockIdx().y == 1
            for val in valid_vals
                @cuprintln(real(val))
            end
        end

        # 4) Get median value
        mid = div(not_nan_count, 2)
        if not_nan_count % 2 == 0
            median = (valid_vals[mid] + valid_vals[mid + 1]) / 2
        else
            median = valid_vals[mid + 1]
        end

        if threadIdx().x == 2 && threadIdx().y == 2 && blockIdx().x == 1 && blockIdx().y == 1
            @cushow real(median)
        end

    end
    return
end

# function bencher(U2, histo, histostd)
function bencher()

    # Test U2 matrix
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