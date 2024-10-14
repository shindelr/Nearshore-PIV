using CUDA
using DelimitedFiles

function localfilt(x::Matrix{Float32}, y::Matrix{Float32}, u::Matrix{Float32}, 
    v::Matrix{Float32}, threshold::Int32, median_bool=true, m=3)

    dim1 = round(Int32, size(u, 1) + 2 * floor(m / 2))
    dim2 = round(Int32, size(u, 2) + 2 * floor(m / 2))

    nu = zeros(eltype(u), (dim1, dim2)) * NaN
    nv = zeros(eltype(u), (dim1, dim2)) * NaN

    # Transfer over data
    from_cols = round(Int32, floor(m / 2) + 1)
    minus_rows = round(Int32, floor(m / 2))
    nu[from_cols:end-minus_rows, from_cols:end-minus_rows] = u
    nv[from_cols:end-minus_rows, from_cols:end-minus_rows] = v

    U2::CuArray{ComplexF32} = nu .+ im .* nv

    ma, na = size(U2)
    histostd = CUDA.zeros(ComplexF32, size(nu))
    histo = CUDA.zeros(ComplexF32, size(nu))

    # Run the kernel

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
