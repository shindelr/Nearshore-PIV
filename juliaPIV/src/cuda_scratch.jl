using CUDA
using DelimitedFiles

function localfilt(x::Matrix{Float32}, y::Matrix{Float32}, u::Matrix{Float32}, 
    v::Matrix{Float32}, threshold::Int32, median_bool=true, m=3)

    dim1 = round(Int32, size(u, 1) + 2 * floor(m / 2))
    dim2 = round(Int32, size(u, 2) + 2 * floor(m / 2))

    # nu = zeros(eltype(u), (dim1, dim2)) * NaN
    # nv = zeros(eltype(u), (dim1, dim2)) * NaN
    nu = CUDA.zeros(Float32, (dim1, dim2)) * NaN
    nv = CUDA.zeros(Float32, (dim1, dim2)) * NaN

    # Transfer over data
    from_cols = round(Int32, floor(m / 2) + 1)
    minus_rows = round(Int32, floor(m / 2))
    nu[from_cols:end-minus_rows, from_cols:end-minus_rows] = u
    nv[from_cols:end-minus_rows, from_cols:end-minus_rows] = v

    INx = CUDA.zeros(Float32, size(nu))

    U2::CuArray{ComplexF32} = nu .+ im .* nv

    ma, na = size(U2)
    histostd = CUDA.zeros(ComplexF32, size(nu))
    histo = CUDA.zeros(ComplexF32, size(nu))

    # Pre-GPU small calculations
    i_bound = na - m + 2
    j_bound = ma - m + 2
    m_floor_two = floor(Int32, m / 2)
    m_ceil_two = ceil(Int32, m / 2)

    # Run the kernel
    # @cuda filt_loop_kernel!(i_bound, j_bound, m_floor_two, m_ceil_two, INx, U2, histo, histostd) 

    for ii in m-1:1:na-m+2
        for jj in m-1:1:ma-m+2
            # Get 3x3 submatrix
            m_floor_two = floor(Int32, m / 2)
            tmp::CuArray{ComplexF32} = U2[jj-m_floor_two:jj+m_floor_two,
                                            ii-m_floor_two:ii+m_floor_two]

            # Set center to NaN for some reason
            tmp[ceil(Int32, m / 2), ceil(Int32, m / 2)] = NaN

            # Run the appropriate stat depending on method arg.
            histo[jj, ii] = median_bool ? im_median_magnitude(tmp[:]) : mean(tmp[:])
            histostd[jj, ii] = im_std(tmp[:])
        end
    end

    # Locate gridpoints w/higher value than the threshold
    coords = findall(
    (real(U2) .> real(histo) .+ threshold .* real(histostd)) .|
    (imag(U2) .> imag(histo) .+ threshold .* imag(histostd)) .|
    (real(U2) .< real(histo) .- threshold .* real(histostd)) .|
    (imag(U2) .< imag(histo) .- threshold .* imag(histostd)))

    # Then "filter" those points out by changing them to NaN!
    for jj in eachindex(coords)
        nu[coords[jj]] = NaN
        nv[coords[jj]] = NaN
    end

    m_ceil_two = ceil(Int32, m / 2)
    m_floor_two = floor(Int32, m / 2)
    hu::Matrix{Float32} = nu[m_ceil_two:end-m_floor_two, m_ceil_two:end-m_floor_two]
    hv::Matrix{Float32} = nv[m_ceil_two:end-m_floor_two, m_ceil_two:end-m_floor_two]

    return hu, hv
end

function filt_loop_kernel!(i_bound, j_bound, m_floor_two, m_ceil_two,
                            INx, U2, histo, histostd)
    # GPUs operate in blocks with many threads in each block.
    # To take advantage of this, we can parallelize the double for loop
    # in the original CPU code.
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    j = threadIdx().y + (blockIdx().y - 1) * blockDim().y

    # Each thread is operating in parallel, so just check that its given
    # indices are within the bound of the matrix.
    if i <= i_bound && j <= j_bound
        if INx[j, 1] != 1
            tmp = CUDA.zeros(Float32, (j+m_floor_two, i+m_floor_two)) * NaN
            tmp = U2[j - m_floor_two: j + m_floor_two,
                      i - m_floor_two: i + m_floor_two]
            tmp[m_ceil_two, m_ceil_two] = NaN
            usum = 1
            # usum = im_median_magnitude(tmp[:])
            # histostd[j, i] = im_std(tmp[:])
        else
            usum = NaN
            tmp = NaN
            histostd[j, i] = NaN
        end
        histo[j, i] = usum
    end
    return
end

function im_median_magnitude(collection::AbstractArray{ComplexF32})
    i = filter(x -> !isnan(x), collection)
    isempty(i) && return NaN
    n = length(i)
    v = partialsort!(i, div(n + 1, 2, RoundDown):div(n + 1, 2, RoundUp); by=x -> (abs2(x), angle(x)))
    return sum(v) / length(v)
end

function im_std(collection::AbstractArray{ComplexF32})
    i = filter(x -> !isnan(x), collection)

    if length(i) > 0
        real_part = std(real(i), corrected=false)
        im_part = std(imag(i), corrected=false)
        return real_part + im_part * im
    end

    # All NaN!
    return NaN
end

function testing()
    println("Running localfilt CUDA test")
    
    x = readdlm("../../tests/gpu_tests/x.csv", ',', Float32)
    y = readdlm("../../tests/gpu_tests/y.csv", ',', Float32)
    u = readdlm("../../tests/gpu_tests/datax.csv", ',', Float32)
    v = readdlm("../../tests/gpu_tests/datay.csv", ',', Float32)
    threshold::Int32 = 3
    
    @show eltype(x) eltype(y) eltype(u) eltype(v)
    datax, datay = localfilt(x, y, u, v, threshold)
end
