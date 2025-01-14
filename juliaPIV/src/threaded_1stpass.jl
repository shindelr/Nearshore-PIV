# include("main.jl")
function threaded_firstpass(A::T, B::T, N::Int32, overlap::Float32,
    idx::Matrix{Float32}, idy::Matrix{Float32}) where {T}

    M = N

    # Set up for FFT plans
    pad_matrix_a = pad_for_xcorr(A[1:M, 1:N])
    # pad_matrix_b = pad_for_xcorr(B[1:M, 1:N])
    P = plan_fft(pad_matrix_a, flags=FFTW.ESTIMATE)
    Pi = plan_ifft(pad_matrix_a, flags=FFTW.ESTIMATE)

    # Initializing matrices
    sy, sx = size(A)
    xx_dim1 = ceil(Int32, ((size(A, 1) - N) / ((1 - overlap) * N))) + 1
    xx_dim2 = ceil(Int32, ((size(A, 2) - M) / ((1 - overlap) * M))) + 1
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))
    yy = zeros(eltype(A), (xx_dim1, xx_dim2))

    # New displacement matrices
    data_dim_1::Int32 = round(Int32, sy / (N * (1 - overlap)))
    data_dim_2::Int32 = round(Int32, sx / (N * (1 - overlap)))
    datax = zeros(Float32, (data_dim_1, data_dim_2))
    datay = zeros(eltype(A), (xx_dim1, xx_dim2))

    num_windows = length(1:((1-overlap)*N):(sy-N+1)) * length(1:((1-overlap)*M):(sx-M+1))
    winds_A = Array{Float32}(undef, N, M, num_windows)
    winds_B = Array{Float32}(undef, N, M, num_windows)
    mf = nextpow(2, N + M)
    pad_a = fill!(Array{ComplexF32}(undef, mf, mf, num_windows), 0.0)
    pad_b = fill!(Array{ComplexF32}(undef, mf, mf, num_windows), 0.0)
    # R = fill!(Array{Float32}(undef, mf, mf, num_windows), 0.0)

    cj = 1
    k = 1
    for jj in 1:((1-overlap)*N):(sy-N+1)  # Iterate in steps of half a window size
        ci = 1 
        for ii in 1:((1-overlap)*M):(sx-M+1) 

            if isnan(idx[cj, ci])
                idx[cj, ci] = 0
            end

            if isnan(idy[cj, ci])
                idy[cj, ci] = 0
            end

            # idx/y are matrices containing pixel displacement data in x/y directions
            # Adjust displacement coordinates for index bounding
            if (jj + idy[cj, ci]) < 1
                idy[cj, ci] = 1 - jj 
            elseif (jj + idy[cj, ci]) > (sy - N + 1)
                idy[cj, ci] = sy - N + 1 - jj 
            end

            if (ii + idx[cj, ci]) < 1
                idx[cj, ci] = 1 - ii
            elseif (ii + idx[cj, ci]) > (sx - M + 1)
                idx[cj, ci] = sx - M + 1 - ii
            end

            # Extract each window from Image A and Image B
            C = A[floor(Int32, jj):floor(Int32, jj + N - 1), 
                floor(Int32, ii):floor(Int32, ii + M - 1)]
            D = B[floor(Int32, jj + idy[cj, ci]):floor(Int32, jj + N - 1 + idy[cj, ci]),
                floor(Int32, ii + idx[cj, ci]):floor(Int32, ii + M - 1 + idx[cj, ci])]

            # Normalize each window against its mean
            winds_A[:, :, k] = (C .- mean(C))
            winds_B[:, :, k] = (D .- mean(D))

            # This bit isn't consequential at the moment
            xx[cj, ci] = ii + M / 2
            yy[cj, ci] = ii + N / 2            

            ci += 1
            k += 1
        end
        cj += 1
    end

    # Cross correlate on the GPU!
    R = multi_xcorr(winds_A, winds_B, P, Pi, pad_a, pad_b)
    
    ci = 1; cj = 1
    for k in axes(R, 3)
        r = R[:, :, k]
        # Find position of maximal value of R
        max_coords = Vector{NTuple{2, Float32}}()
        subset = r[Int32(0.5 * N + 2):Int32(1.5 * N - 3), Int32(0.5 * M + 2):Int32(1.5 * M - 3)]
        fast_max!(max_coords, subset)

        # Adjust for subset positions
        max_coords = [(i[1] + Int32(0.5 * N + 1), 
                        i[2] + Int32(0.5 * M + 1))
                        for i in max_coords]

        # Handle a vector that has multiple maximum coordinates. Take the 
        # weighted average of the coordinates.
        if length(max_coords) > 1
            max_x1 = round(Int32, sum([c[2] * i for (i, c) in enumerate(max_coords)]) / sum([c[2] for c in max_coords]))
            max_y1 = round(Int32, sum([c[1] * i for (i, c) in enumerate(max_coords)]) / sum([c[1] for c in max_coords]))

        elseif isempty(max_coords)
            idx[cj, ci] = NaN
            idy[cj, ci] = NaN
            max_x1 = NaN
            max_y1 = NaN

        # Otherwise, unpack into max coordinates
        else
            max_y1, max_x1 = max_coords[1][1], max_coords[1][2]
        end

        # Store displacements in variables datax/datay
        datax[cj, ci] -= (max_x1 - M) + idx[cj, ci]
        datay[cj, ci] -= (max_y1 - M) + idy[cj, ci]
        xx[cj, ci] = k + M / 2
        yy[cj, ci] = k + N / 2

        ci += 1
        if ci > size(datax, 2)
            ci = 1
            cj += 1
        end
    end

return xx, yy, datax, datay
end

function multi_xcorr(A::Array{Float32}, B::Array{Float32}, 
    plan::FFTW.cFFTWPlan, iplan::AbstractFFTs.ScaledPlan,
    pad_matrix_a::Array{ComplexF32}, pad_matrix_b::Array{ComplexF32})

    N = size(A, 3)
    R = Array{Float32}(undef, size(pad_matrix_a, 1) - 1, size(pad_matrix_a, 2) - 1, N)

    Threads.@threads for k in 1:N
        R[:, :, k] = xcorrf2(A[:, :, k], B[:, :, k], plan, iplan, 
                        pad_matrix_a[:, :, k], pad_matrix_b[:, :, k])
    end

    return R
end
