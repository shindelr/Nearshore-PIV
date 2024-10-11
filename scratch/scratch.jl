using CUDA
using StaticArrays
using DelimitedFiles

function localfilt_kernel!(U2, histo, histostd)
    i = threadIdx().x + (blockIdx().x - Int32(1)) * blockDim().x
    j = threadIdx().y + (blockIdx().y - Int32(1)) * blockDim().y


    if i > 1 && i < size(U2, 1) && j > 1 && j < size(U2, 2)
        # 1) Get subarray info
        @inbounds begin
            subarray = @view U2[1:3, 1:3]

            # 2) Filter out NaN values
            valid_vals = @MVector fill(ComplexF32(NaN + NaN * im), 9)

            not_nan_count = 0
            for val in subarray
                if !isnan(real(val)) && !isnan(imag(val)) 
                    not_nan_count += 1
                    valid_vals[not_nan_count] = val
                end
            end
            if not_nan_count == 0
                # histo[j, i] = NaN Might not be any need because it's already NaN
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

            # 4) Get median value
            mid = div(not_nan_count, 2)
            if not_nan_count % 2 == 0
                median = (valid_vals[mid] + valid_vals[mid + 1]) / 2
            else
                median = valid_vals[mid + 1]
            end

            # 5) Get standard deviation
            real_mean = 0
            im_mean = 0
            real_var = 0
            im_var = 0
            for k in 1:not_nan_count
                real_mean += real(valid_vals[k])
                im_mean += imag(valid_vals[k])
            end
            real_mean /= not_nan_count
            im_mean /= not_nan_count

            for k in 1:not_nan_count
                real_var += (valid_vals[k] - real_mean)^2
                im_var += (valid_vals[k] - real_mean)^2
            end
            real_std = sqrt(real_var / not_nan_count)
            im_std = sqrt(im_var / not_nan_count)

            std = real_std + im_std * im

            # 6) Set histo[j, i] and histostd[j, i]
            histo[i, j] = median
            histostd[i, j] = std
        end

    end
    return
end

# function localfilt_kernel_config(U2)
function localfilt_kernel_config()
    U2 = CUDA.CuArray(ComplexF32[
        3 + 4im 1 + 2im 7 + 8im 0 + 0im 5 + 6im;
        9 + 10im 11 - 12im 15 - 16im 13 + 14im 7 - 8im;
        21 + 22im 17 + 18im 19 - 20im 25 + 26im 15 - 16im;
        27 - 28im 23 - 24im 29 + 30im 31 - 32im 25 + 26im
    ])

    histo = CUDA.zeros(ComplexF32, size(U2, 1), size(U2, 2))
    histostd = CUDA.zeros(ComplexF32, size(U2, 1), size(U2, 2))

    # 2D array size
    M = size(U2, 1)
    N = size(U2, 2)

    # Get config for kernel
    kernel = @cuda launch=false localfilt_kernel!(U2, histo, histostd)
    config = launch_configuration(kernel.fun)

    # Specificy nukber of threads and blocks for 2D array
    threads_x = Int32(min(M, sqrt(config.threads)))
    threads_y = Int32(min(N, sqrt(config.threads)))
    # For aquilla, the threads cannot be more than 1024
    @assert threads_x * threads_y <= config.threads

    blocks_x = Int32(cld(M, threads_x))
    blocks_y = Int32(cld(N, threads_y))
    
    CUDA.@sync begin
        kernel(U2, histo, histostd; 
               threads=(threads_x, threads_y), blocks=(blocks_x, blocks_y)
            )
    end
    return histo
end

function localfilt(U2)
    # U2 = @SMatrix ComplexF32[
    #     3 + 4im 1 + 2im 7 + 8im 0 + 0im 5 + 6im;
    #     9 + 10im 11 - 12im 15 - 16im 13 + 14im 7 - 8im;
    #     21 + 22im 17 + 18im 19 - 20im 25 + 26im 15 - 16im;
    #     27 - 28im 23 - 24im 29 + 30im 31 - 32im 25 + 26im
    # ]
 
    histo = zeros(ComplexF32, size(U2, 1), size(U2, 2))
    histostd = zeros(ComplexF32, size(U2, 1), size(U2, 2))

    ma, na = size(U2)
    m = 3

    for ii in m-1:1:na-m+2
        for jj in m-1:1:ma-m+2
            # Get a 3x3 submatrix of U2
            m_floor_two = floor(Int32, m / 2)
            tmp = U2[jj-m_floor_two:jj+m_floor_two,
                        ii-m_floor_two:ii+m_floor_two]

            # Assign the center value to NaN
            tmp[ceil(Int32, m / 2), ceil(Int32, m / 2)] = NaN

            # Run the appropriate stat depending on method arg.
            histo[jj, ii] = im_median_magnitude(tmp[:])
            histostd[jj, ii] = im_std(tmp[:])
        end
    end
    return histo, histostd
end

function im_median_magnitude(collection)
    i = filter(x -> !isnan(x), collection)
    isempty(i) && return NaN
    n = length(i)
    v = partialsort!(i, 
                    div(n + 1, 2, RoundDown):div(n + 1, 2, RoundUp); 
                    by=x -> (abs2(x), angle(x))
                )
    return sum(v) / length(v)
end

function im_std(collection)
    i = filter(x -> !isnan(x), collection)

    if length(i) > 0
        real_part = std(real(i), corrected=false)
        im_part = std(imag(i), corrected=false)
        return real_part + im_part * im
    end

    # All NaN!
    return NaN
end

function build_cu_u2()
    nv = readdlm("../tests/gpu_tests/nv.csv", ',')
    nu = readdlm("../tests/gpu_tests/nu.csv", ',')

    nv_cu = CuArray{Float32}(undef, size(nv))
    nu_cu = CuArray{Float32}(undef, size(nu))

    copyto!(nv_cu, nv)
    copyto!(nu_cu, nu)

    U2 = nu_cu .+ im .* nv_cu
    localfilt_kernel_config(U2)
end

function build_normal_u2()
    nv = readdlm("../tests/gpu_tests/nv.csv", ',')
    nu = readdlm("../tests/gpu_tests/nu.csv", ',')

    U2::Matrix{ComplexF32} = nu .+ im .* nv
    histo, histostd = localfilt(U2)
end