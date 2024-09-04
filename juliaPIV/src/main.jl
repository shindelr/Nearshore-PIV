# Third party modules
using Statistics
using FFTW            # Fast Fourier Transforms library built on C
using Images          # Basic image processing library
using FileIO          # I/O library
using DelimitedFiles  # Write matrices to CSV
using Interpolations
using Plots
using Luxor            # For creating inpolygon() functionality

# PASS FUNCTIONS 
"""
### multipassx
    \n**:params:**\n
    A: Matrix containing image data of first frame.\n
    B: Matrix containing image data of second frame.\n
    wins: 2D matrix of ints containing sub-window pixel sizes for each pass.\n
    Dt: Frame time step in seconds (int). Pass 1 for 'pixels per frame velocity'.\n
    overlap: Fraction of window overlap. Int.\n
    sensit: Threshold for vector validation. Int. \n
    \n**:returns:**\n
    x: \n
    y: \n
    u: \n
    v: \n
    SnR: Ratio representing signal-to-noise.\n
    Pkh: Peak height for use in validation of vector field?\n
"""
function multipassx(A::Matrix{Float32}, B::Matrix{Float32}, wins::Vector{Int32}, Dt::Int32, 
                    overlap::Float32, sensit::Int32)
    sy, sx = size(A)
    # total_passes = size(wins, 1)
    total_passes = length(wins)

    # Initial passes are for removing large-scale displacements. Initialize
    # displacements (datax,datay) to zero
    data_dim_1::Int32 = sy/(wins[1] * (1-overlap))
    data_dim_2::Int32 = sx/(wins[1] * (1-overlap))
    datax = zeros(Float32, (data_dim_1, data_dim_2))
    datay = copy(datax)

    for i in 1:total_passes - 1
        println("Pass ", i, " of ", total_passes )
    
        x, y, datax, datay = firstpass(A, B, wins[i], overlap, datax, datay)
        datax, datay = localfilt(x, y, datax, datay, sensit)
        
        datax, datay = linear_naninterp(datax, datay)
        datax = floor.(datax)
        datay = floor.(datay)

        if i != total_passes - 1
            Y, X, YI, XI = build_grids_2(datax)
            datax = regular_interp(datax, X, Y, XI, YI)
            datay = regular_interp(datay, X, Y, XI, YI)
            
            # NAN BORDER ATTEMPT ----------------
            datax = make_nan_border(datax)
            datay = make_nan_border(datay)

            datax, datay = linear_naninterp(datax, datay)
            datax = round.(datax)
            datay = round.(datay)
            # ---------------------------------
        end
    end

    println("Final Pass")
    x, y, u, v, SnR, Pkh = finalpass(A, B, wins[end], overlap, datax, datay, Dt)

    return x, y, u, v, SnR, Pkh
end

"""
### firstpass
    *Note: First pass is a misnomer for this function, as it's called N-1 times
    depending on the size of the given windows. Consider renaming*\n\n
    Set up matrix indices along the image frame according to the desired overlap and \
    window sizes. Calls **xcorrf2** which finally uses FFTs to calculate cross \
    correlation.
    \n**:params:**\n
    A: Matrix containing image data of first frame.\n
    B: Matrix containing image data of second frame.\n
    N: Pair of Float64 representing sub-window pixel sizes for the pass.\n
    overlap: Fraction of window overlap. Int.\n
    idx: Matrix of same type as A, containing data displacement information.\n
    idy: Matrix of same type as A, containing data displacement information.\n
    \n**:returns:**\n
    x: \n
    y: \n
    datax: \n
    datay: \n
"""
function firstpass(A::Matrix{Float32}, B::Matrix{Float32}, N::Int32, overlap::Float32, 
                    idx::Matrix{Float32}, idy::Matrix{Float32})
    M = N

    # Set up for FFT plans
    pad_matrix_a = pad_for_xcorr(A[1:M, 1:N])
    pad_matrix_b = pad_for_xcorr(B[1:M, 1:N])
    P = plan_fft(pad_matrix_a; flags=FFTW.PATIENT)
    Pi = plan_ifft(pad_matrix_a; flags=FFTW.PATIENT)

    # Initializing matrices
    sy, sx = size(A)
    xx_dim1 = ceil(Int32, ((size(A,1)-N) / ((1-overlap) * N))) + 1
    xx_dim2 = ceil(Int32, ((size(A,2)-M) / ((1-overlap) * M))) + 1
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))
    yy = copy(xx)
    datax = copy(xx)
    datay = copy(xx)
    IN = zeros(Int32, size(A))

    cj = 1

    for jj in 1:((1-overlap) * N):(sy - N + 1)
        ci = 1
        for ii in 1:((1-overlap) * M):(sx - M + 1)
            # Floor correct?
            IN_i_1 = floor(Int32, (jj + N/2))
            IN_i_2 = floor(Int32, (ii + M/2))

            if IN[IN_i_1, IN_i_2] != 1

                if isnan(idx[cj, ci])
                    idx[cj, ci] = 0
                end

                if isnan(idy[cj, ci])
                    idy[cj, ci] = 0
                end

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

                C = A[floor(Int32, jj):floor(Int32, jj+N-1), 
                      floor(Int32, ii):floor(Int32, ii+M-1)]
                D = B[floor(Int32, jj+idy[cj, ci]):floor(Int32, jj+N-1+idy[cj, ci]), 
                      floor(Int32, ii+idx[cj, ci]):floor(Int32, ii+M-1+idx[cj, ci])]
                
                C = C.-mean(C); D = D.-mean(D)
                stad1 = std(C); stad2 = std(D)

                if stad1 == 0
                    stad1 = NaN
                end
                if stad2 == 0
                    stad2 = NaN
                end

                # Call xcorrf2, passing in the FFT plan and normalize result
                R = xcorrf2(C, D, P, Pi, pad_matrix_a, pad_matrix_b) ./ ( N * M * stad1 * stad2)

                # Find position of maximal value of R
                if size(R, 1) == (N - 1)  # I think checking for second to last pass
                    max = maximum(R)
                    max_coords = findall(x -> x == max, R)
                else
                    subset = R[Int32(.5*N+2):Int32(1.5*N-3), 
                                Int32(.5*M+2):Int32(1.5*M-3)]
                    max = maximum(subset)
                    max_coords = findall(x -> x == max, subset)
                    max_coords = [(i[1] + Int32(0.5*N+1), i[2] + Int32(0.5*M+1)) for i in max_coords]
                end

                # Handle a vector that has multiple maximum coordinates.
                # Sum the product of each x and y indice with its own indice within
                # the max_coords vector.
                if length(max_coords) > 1
                    # FLIPPED C,I TO I,C
                    max_x1 = round(
                        Int32, sum([c[2] * i for (i, c) in enumerate(max_coords)]) / 
                        sum([c[2] for c in max_coords]))
                    max_y1 = round(
                        Int32, sum([c[1] * i for (i, c) in enumerate(max_coords)]) / 
                        sum([c[1] for c in max_coords]))
                    
                # Handle empty max_coords vector.
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
                datax[cj, ci] -= (max_x1 - M) + idx[cj,ci]
                datay[cj, ci] -= (max_y1 - M) + idy[cj,ci]
                xx[cj, ci] = ii + M/2
                yy[cj, ci] = jj + N/2
                ci += 1

            else
                xx[cj, ci] = ii + M/2
                yy[cj, ci] = jj + N/2
                datax[cj, ci] = NaN
                datay[cj, ci] = NaN
                ci += 1
            end
        end

        cj += 1
    end
    return xx, yy, datax, datay
end

"""
### finalpass
    TODO: Write me.

    OG docstring
    Provides the final pass to get the displacements with
    subpixel resolution.



    1999 - 2011, J. Kristian Sveen (jks@math.uio.no)
    For use with MatPIV 1.7, Copyright
    Distributed under the terms of the GNU - GPL license
    timestamp: 09:26, 4 Mar 2011
"""
function finalpass(A::Matrix{Float32}, B::Matrix{Float32}, N::Int32, ol::Float32, 
                    idx::Matrix{Float32}, idy::Matrix{Float32}, Dt::Int32)
    M = N

    # FFT setup
    pad_matrix_a = pad_for_xcorr(A[1:M, 1:N])
    pad_matrix_b = pad_for_xcorr(B[1:M, 1:N])
    P = plan_fft(pad_matrix_a; flags=FFTW.PATIENT)
    Pi = plan_ifft(pad_matrix_a; flags=FFTW.PATIENT)

    cj = 1
    sy, sx = size(A)
    dim_1 = ceil(Int32, (sy - N) / ((1 - ol) * N)) + 1
    dim_2 = ceil(Int32, (sx - M) / ((1 - ol) * M)) + 1

    xp = zeros(eltype(A), (dim_1, dim_2))
    yp = copy(xp); up = copy(xp); vp = copy(xp); SnR = copy(xp); Pkh = copy(xp)
   
    # Main pass loop
    for jj in 1:((1 - ol) * N):sy - N + 1
        ci = 1
        for ii in 1:((1 - ol) * M):sx - M + 1

            if isnan(idx[cj, ci])
                idx[cj, ci] = 0
            end

            if isnan(idy[cj, ci])
                idy[cj, ci] = 0
            end

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

            D2 = B[
                Int32(jj + idy[cj, ci]):Int32(jj + N - 1 + idy[cj, ci]),
                Int32(ii + idx[cj, ci]):Int32(ii + M - 1 + idx[cj, ci])
                ]
            E = A[Int32(jj):Int32(jj + N - 1),
                  Int32(ii):Int32(ii + M - 1)]

            stad1 = std(E); stad2 = std(D2)

            if stad1 == 0
                stad1 = 1
            end
            if stad2 == 0
                stad2 = 1
            end

            E = E.-mean(E); F = D2.-mean(D2)

            # Cross correlate and FFT
            R = xcorrf2(E, F, P, Pi, pad_matrix_a, pad_matrix_b) ./ ( N * M * stad1 * stad2)

            if !any(isnan.(R)) & !all(x -> x == 0, R)

                # Find position of maximal value of R
                if size(R, 1) == (N - 1)
                    max_coords = findall(x -> x == maximum(R), R)
                else
                    subset = R[Int32(.5 * N + 2):Int32(1.5 * N - 3), 
                                Int32(.5 * M + 2):Int32(1.5 * M - 3)]

                    # max_coords::Vector{Tuple{Int32, Int32}} = findall(x -> x == maximum(subset), subset)
                    max_coords = findall(x -> x == maximum(subset), subset)
                    max_coords = [(i[1] + Int32(0.5*N+1), i[2] + Int32(0.5*M+1)) for i in max_coords]
                end

                # Handle a vector that has multiple maximum coordinates.
                # Sum the product of each x and y indice with its own indice within
                # the max_coords vector.
                if length(max_coords) > 1
                    @show max_coords
                    max_x1 = round(Int32, 
                            sum([c[2]^2 for c in max_coords]) /
                            sum([c[2] for c in max_coords]))
                    max_y1 = round(Int32, 
                            sum([c[1]^2 for c in max_coords]) /
                            sum([c[1] for c in max_coords]))
                end
                
                # Unpack cartesian index type. Only a handful iterations here
                if length(max_coords) == 1
                    max_y1, max_x1 = max_coords[1][1], max_coords[1][2]
                end
                
                # Some kind of manual adjustment?
                if max_x1 == 1
                    max_x1 = 2
                end
                if max_y1 == 1
                    max_y1 = 2
                end

                # Runs without error, seems to be equivalent.
                # 3-point peak fit using gaussian fit
                x_0, y_0 = intpeak(max_x1, max_y1, 
                                R[max_y1, max_x1],
                                R[max_y1, max_x1 - 1],
                                R[max_y1, max_x1 + 1],
                                R[max_y1 - 1, max_x1],
                                R[max_y1 + 1, max_x1],
                                N
                )

                R2 = copy(R)

                # This section had a note to try to simplify their try-catch
                # clause by using a distance check.
                try
                    R2[max_y1 - 3:max_y1 + 3, max_x1 - 3: max_x1 + 3] .= NaN
                catch
                    R2[max_y1 - 1: max_y1 + 1, max_x1 - 1: max_x1 + 1] .= NaN
                end

                if size(R, 1) == (N - 1)
                    max_val = maximum(R2)
                    p2_coords = findall(x -> x == max_val, R2)
                    if length(p2_coords) == 1
                        p2_y2, p2_x2 = p2_coords[1][1], p2_coords[1][2]
                    end
                else
                    # Find subset of R2 and remember where all !NaN vals are for later reindexing
                    subset = R2[Int32(0.5 * N):Int32(1.5 * N - 1), 
                                Int32(0.5 * M):Int32(1.5 * M - 1)]
                    non_nans = findall(x -> !isnan(x), subset)

                    # Filter out all non-NaN values to find the maximum
                    # NaN vals contaminate arrays for the min/max functions
                    filt_subset = filter(!isnan, subset)
                    subset_max_val = maximum(filt_subset)
                    max_coords_2 = findall(x -> x == subset_max_val, filt_subset)
                    
                    # Re-adjust coordinates to match original coords in R2
                    og_max_coords = non_nans[max_coords_2]
                    max_coords_2 = [(i[1] + floor(Int32, 0.5 * N - 1), 
                                   i[2] + floor(Int32, 0.5 * M - 1))
                                   for i in og_max_coords]
                end

                if length(max_coords_2) == 1
                    p2_y2, p2_x2 = max_coords_2[1][1], max_coords_2[1][2]
                elseif length(max_coords_2) > 1
                    p2_x2 = round(Int32, length(max_coords_2) ./ 2)
                    p2_y2 = round(Int32, length(max_coords_2) ./ 2)
                elseif isempty(max_coords)
                    error("Empty set found in Final Pass")
                end

                snr = R[max_y1, max_x1] / R2[p2_y2, p2_x2]
                SnR[cj, ci] = snr
                up[cj, ci] = (-x_0 + idx[cj, ci]) / Dt
                vp[cj, ci] = (-y_0 + idy[cj, ci]) / Dt
                xp[cj, ci] = ii + (M / 2) - 1
                yp[cj, ci] = jj + (N / 2) - 1
                Pkh[cj, ci] = R[max_y1, max_x1]
                
            else
                up[cj, ci] = NaN
                vp[cj, ci] = NaN
                SnR[cj, ci] = NaN
                Pkh[cj, ci] = 0
                xp[cj, ci] = ii + M / 2 - 1
                yp[cj, ci] = jj + N / 2 - 1
            end
            ci += 1
        end
        cj += 1
    end
    return xp, yp, up, vp, SnR, Pkh
end


# FOURIER
"""
### pad_for_xcorr
    Determine the dimensions of a padded matrix to be used for xcorrf2. Depends
    on the dimensions of the current window size and the original array. Takes
    the current window size, scales it by a power of 2, then creates an array of 
    zeros that size. 
    :params:
        - trunc_matrix: The original matrix to be padded, but truncated down to\
        the window's size. \n
    :returns:
        - A padded matrix of zeros.
"""
function pad_for_xcorr(trunc_matrix::Matrix{Float32})
    ma, na = size(trunc_matrix)
    mf = nextpow(2, ma + na)  
    return zeros(eltype(trunc_matrix), (mf, mf)) 
end

"""
### xcorrf2
    Two-dimensional cross-correlation using Fourier transforms.
    XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
    If desired, you can pass in `None` for the parameter `plan`, opting to\
    us the original "padding" method provided in MatPIV. This method is not \
    nearly as effective in Julia and it's recommended that you use the plan.
    \n**:params:**\n
    A: matrix (2D array) to be compared.\n
    B: matrix ((2D array)) to be compared.\n
    plan: Callable function representing a pre-planned, omptimized FFT. Created \
    using the dimensions of matrices A & B.
    \n**:return:**\n
    c: A matrix whose values reflect the 2D correlation between every cell \
    in A & B. Matrix is 2D float 64.

    Originally written in Matlab by,\n
    Author(s): R. Johnson\n
    Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(A, B, plan, iplan, pad_matrix_a, pad_matrix_b)
    # Unpack size() return tuple into appropriate variables
    ma, na = size(A)
    mb, nb = size(B)
    
    # Reverse conjugate
    B = conj(B[mb:-1:1, nb:-1:1])

    # Transfer data from og matrix to optimized sized ones
    pad_matrix_a[1:size(A,1), 1:size(A,2)] = A[1:size(A,1), 1:size(A,2)]
    pad_matrix_b[1:size(B,1), 1:size(B,2)] = B[1:size(B,1), 1:size(B,2)]

    # Runs optimized FFTs using the plan
    at = plan * pad_matrix_b
    bt = plan * pad_matrix_a

    # Mult transforms and invert
    mult_at_bt = at.*bt
    c = iplan * mult_at_bt
    # c = ifft(mult_at_bt)

    # Make all real
    c = real(c)

    # Trim
    rows = ma + mb
    cols = na + nb
    # Keep everything from the first index to the index where padding began
    c = c[1:rows - 1, 1:cols - 1]
    return c
end


# FILTERS and Interpolations
"""
### linear_naninterp
    Interpolates NaN's in a vectorfield. Sorts all spurious 
    vectors based on the number of spurious neighbors to a 
    point. The function replaces NaN values in the input 
    vector fields U and V using linear interpolation.\n
    Interpolation starts with the ones that have the least 
    number of outliers in their neighborhood and loops until no 
    NaN's are present in the field.\n
    NOTE: This function is a completely gutted version of the
    original naninterp combined with naninterp2. The majority
    of naninterp was varargin code figuring out which method
    to execute naninterp2 with. Since arguments to naninterp
    were hardcoded in, the desired form of naninterp was just
    turned into a specific function. If necessary, other 
    methods can be written to emulate the functionality of
    the original functions. Moreover, though a mask was being
    passed in, the masking branch was not being executed in
    either functions. Therefore, in this implementation, masking
    was left out.\n

    Parameters:
    -----------
    - u, v: `Matrices`\n

    Original author:
    ----------------
    J. Kristian Sveen (jks@math.uio.no)
    Department of Mathematics, Mechanics Division, University of Oslo, Norway
    Copyright 1999 - 2001
    For use with MatPIV 1.6, Copyright
    Distributed under the terms of the GNU - GPL license
"""
function linear_naninterp(u::Matrix{Float32}, v::Matrix{Float32})
    coords::Vector{CartesianIndex{2}} = findall(x->isnan(x), u)
    numm = length(coords)
    dy,dx = size(u)
    lp = 1; tel = 1

    # Now sort the NaN's after how many neighbors they have that are
    # physical values. Then we first interpolate those that have 8
    # neighbors, followed by 7, 6, 5, 4, 3, 2 and 1
    while !isempty(coords)
        nei = zeros(Int64, length(coords), 3)

        # Check neighbors
        for i in eachindex(coords)
            py = coords[i][1]; px = coords[i][2]
            corx1 = 0; corx2 = 0; cory1 = 0; cory2 = 0

            # Correct if vector is on edge of matrix
            # These are edge cases as we explore each NaN in u.
            if py == 1
                cory1 = 1; cory2 = 0
            elseif py == dy 
                cory1 = 0; cory2 = -1
            end
            if px == 1
                corx1 = 1; corx2 = 0
            elseif px == dx
                corx1 = 1; corx2 = -1
            end
            
            # Create a matrix of NaN's 8 neighbors
            ma = u[
                py - 1 + cory1: py + 1 + cory2,
                px - 1 + corx1: px + 1 + corx2
            ]

            nei[i, 1] = count(!isnan, ma)
            nei[i, 2] = px
            nei[i, 3] = py
        end

        # Sort NEI by row to interpolate vectors with fewest spurious neighbors.
        nei = sortslices(nei, dims=1, lt=Base.isgreater)
  
        # Reconstruct sorted outliers and interpolate 1st 50%.
        idx = findall(x -> x >= 8, nei[:, 1])
        while isempty(idx)
            idx = findall(x -> x >= (8-tel), nei[:, 1])
            tel += 1
        end
        tel = 1
        py = nei[idx, 3]
        px = nei[idx, 2]

        for j in axes(py, 1)
            corx1=0; corx2=0; cory1=0; cory2=0
            if py[j] == 1
                cory1 = 1; cory2 = 0
            elseif py[j] == dy
                cory1 = 0; cory2 = -1
            end
            if px[j] == 1
                corx1 = 1; corx2 = 0
            elseif px[j] == dx
                corx1 = 0; corx2 = -1
            end

            temp_u::Matrix{Float32} = u[py[j] - 1 + cory1:py[j] + 1 + cory2, 
                    px[j] - 1 + corx1:px[j] + 1 + corx2] 
            temp_v::Matrix{Float32} = v[py[j] - 1 + cory1:py[j] + 1 + cory2, 
                    px[j] - 1 + corx1:px[j] + 1 + corx2] 
            
            mean_prep_u = filter(!isnan, temp_u[:])
            mean_prep_v = filter(!isnan, temp_v[:])
            u[py[j], px[j]] = mean(mean_prep_u)
            v[py[j], px[j]] = mean(mean_prep_v)
            
            if lp > numm
                u[py[j], px[j]] = 0
                v[py[j], px[j]] = 0
            end
        end
        coords = findall(x->isnan(x), u)  # Might not be necessary
        lp += 1
    end
    return u, v
end

# Not currently used/functioning
function naninterp(sample, pass)
    nan_coords = findall(x -> isnan(x), sample)
    non_nan_coords= findall(x -> !isnan(x), sample)
    non_nan_coords_matrix = hcat([i[1] for i in non_nan_coords], [i[2] for i in non_nan_coords])
    non_nan_vals= [sample[c] for c in non_nan_coords]
    println("made it: 479") # Last place it gets past

    # Debugging
    if pass == 2
        # writedlm("tests/juliaOut/second_pass_interp/stuck_sample.csv", sample, ',')
        display(non_nan_coords_matrix)
        display(non_nan_vals)
    end

    itp = ScatteredInterpolation.interpolate(InverseMultiquadratic(), non_nan_coords_matrix', non_nan_vals)
    for c in nan_coords
        c_extracted = [c[1]; c[2]]
        itp_val_vec = ScatteredInterpolation.evaluate(itp, c_extracted)
        itp_val = itp_val_vec[1]
        sample[c] = itp_val
    end

    return sample
end

"""
        regular_interp(samples, xs, ys, XI, YI)

    Interpolates the given `samples` using linear interpolation on a regular grid defined by `xs` and `ys`.
    The interpolated values are evaluated at the points defined by `XI` and `YI`.

    Arguments
    ---------
    - `samples::AbstractArray`: The array of samples to be interpolated.
    - `xs::AbstractArray`: The array of x-coordinates defining the grid.
    - `ys::AbstractArray`: The array of y-coordinates defining the grid.
    - `XI::AbstractArray`: The array of x-coordinates where the interpolated values are evaluated.
    - `YI::AbstractArray`: The array of y-coordinates where the interpolated values are evaluated.

    Returns
    ---------
    - `itp_results::Array{Float64, 2}`: The interpolated values evaluated at the points defined by `XI` and `YI`.
"""
function regular_interp(samples, xs, ys, XI, YI)
    itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
    itp_results = [itp(yi, xi) for yi in YI, xi in XI]
    return itp_results
end

"""
        build_grids_2(data)

    Builds coarse and fine grids based on the given data.

    Arguments
    ---------
    `data`: A 2-dimensional array representing the data.

    Returns
    ---------
    `coarse_ys`: A 1-dimensional array representing the coarse grid in the y-direction.
    `coarse_xs`: A 1-dimensional array representing the coarse grid in the x-direction.
    `fine_YI`: A 1-dimensional array representing the fine grid in the y-direction.
    `fine_XI`: A 1-dimensional array representing the fine grid in the x-direction.
"""
function build_grids_2(data)
    coarse_y_dim = size(data, 1)
    coarse_x_dim = size(data, 2)

    min_y, max_y = minimum(data[:, 1]), maximum(data[:, 1])
    min_x, max_x = minimum(data[1, :]), maximum(data[1, :])
    coarse_ys = LinRange(min_y, max_y, coarse_y_dim)
    coarse_xs = LinRange(min_x, max_x, coarse_x_dim)

    # OG no NaN border
    # fine_yi_dim = (coarse_y_dim * 2) + 1
    # fine_xi_dim = (coarse_x_dim * 2) + 1

    # NaN border adjustment
    fine_yi_dim = (coarse_y_dim * 2) - 1
    fine_xi_dim = (coarse_x_dim * 2) - 1

    fine_YI = LinRange(min_y, max_y, fine_yi_dim)
    fine_XI = LinRange(min_x, max_x, fine_xi_dim)

    return coarse_ys, coarse_xs, fine_YI, fine_XI
end

"""
        build_grids(wins, overlap, sx, sy, i)

    Build grids for interpolation.

    Arguments
    ---------
    - `wins::Array{Float64,2}`: A matrix containing the window sizes for each iteration.
    - `overlap::Float64`: The overlap between adjacent windows.
    - `sx::Int`: The size of the x-axis.
    - `sy::Int`: The size of the y-axis.
    - `i::Int`: The index of the current window.

    Returns
    ---------
    - `X::Array{Float64,1}`: The x-coordinates of the coarse grid points.
    - `Y::Array{Float64,1}`: The y-coordinates of the coarse grid points.
    - `XI::Array{Float64,1}`: The x-coordinates of the fine interpolated grid points.
    - `YI::Array{Float64,1}`: The y-coordinates of the fine interpolated grid points.
"""
function build_grids(wins, overlap, sx, sy, i)
        next_win_x = wins[i + 1, 1]
        next_win_y = wins[i + 1, 2]

        # Final window size is duplicated, so check for equality.
        if wins[i, 1] != next_win_x
            X = (1:((1 - overlap) * 2 * next_win_x):
                    sx - 2 * next_win_x + 1) .+ next_win_x
            XI = (1:((1 - overlap) * next_win_x):
                    sx - next_win_x + 1) .+ (next_win_x / 2)
        else
            X = (1:((1 - overlap) * next_win_x):
                    sx - next_win_x + 1) .+ (next_win_x / 2)
            XI = (1:((1 - overlap) * next_win_x):
                    sx - next_win_x + 1) .+ (next_win_x / 2)
            X = copy(XI)
        end

        if wins[i, 2] != next_win_y
            Y = (1:((1 - overlap) * 2 * next_win_y): 
                    sy - 2 * next_win_y + 1) .+ next_win_y
            YI = (1:((1 - overlap) * next_win_y):
                    sy - next_win_y + 1) .+ (next_win_y / 2)
        else
            Y = (1:((1 - overlap) * next_win_y):
                    sy - next_win_y + 1) .+ (next_win_y / 2)
            YI = (1:((1 - overlap) * next_win_y):
                    sy - next_win_y + 1) .+ (next_win_y / 2)
            Y = copy(YI)
        end

    return X, Y, XI, YI
end

function make_border(data)
    # Allocate space for the bordered matrix
    bordered_matrix = zeros(Float64, (size(data, 1) + 2, size(data, 2) + 2))

    # Extract the first and last rows and columns
    first_row = data[1, :]
    last_row = data[end, :]
    first_col = data[:, 1]
    last_col = data[:, end]

    # Place border around new matrix, keeping space for each upcoming new row
    # and column
    bordered_matrix[1, 2:end-1] = first_row
    bordered_matrix[end, 2:end-1] = last_row
    bordered_matrix[2:end-1, 1] = first_col
    bordered_matrix[2:end-1, end] = last_col

    # Fill in the interior of the matrix
    bordered_matrix[2:end-1, 2:end-1] = data
    return bordered_matrix
end

function make_nan_border(data)
    # Allocate space for the bordered matrix, make it NaNs then replace them
    bordered_matrix::Matrix{Float32} = fill(NaN, (size(data, 1) + 2, size(data, 2) + 2))

    # Fill in the interior of the matrix
    bordered_matrix[2:end-1, 2:end-1] = data
    return bordered_matrix
end

"""
        globfilt(u, v)

    Global histogram operator. Find the maximum and minimum velocities allowed 
    vector fields u and v. 

    Arguments
    ---------
        `u`: A 2-dimensional array representing the u-component of the vector field.
        `v`: A 2-dimensional array representing the v-component of the vector field.

    Returns
    ---------
        `u`: A 2-dimensional array representing the filtered u-component.
        `v`: A 2-dimensional array representing the filtered v-component.

    Original Author:
    1999 -2014 copyright J.K.Sveen jks@math.uio.no
    Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
"""
function globfilt(u::Matrix{Float32}, v::Matrix{Float32})
    println("Global filter running - with limit: 3 * std [U V]")
    # Opportunity for perfomance improvement here.
    nan_filt_u = filter(!isnan, u)
    nan_filt_v = filter(!isnan, v)
    scale_test = maximum(sqrt.(nan_filt_u[:].^2 + nan_filt_v[:].^2))
    if scale_test != 0
        scale = 2 / scale_test
    else
        scale = 0.1
    end

    sx = 3 * nan_std(u[:])
    sy = 3 * nan_std(v[:])
    xo = nan_mean(u[:])
    yo = nan_mean(v[:])
    ii = [xo + sx; xo + sx; xo - sx; xo - sx]
    jj = [yo + sy; yo - sy; yo - sy; yo + sy]

    # Locate points inside chosen area
    poly_points = Point.(ii, jj)  # Define polygon
    points = vec(Point.(u, v))    # Convert u, v to Points for isinside
    in = [isinside(p, poly_points) for p in points]  # Check if points are inside polygon

    u[.!in] .= NaN
    v[.!in] .= NaN

    return u, v
end

"""
### localfilt
    Filter out vectors that deviate from the median or the mean of their \
    surrounding neighbors by the factor `threshold` times the standard \
    deviation of the neighbors.\n
    Parameters:
    -----------
    - x, y, u, v : `Matrices`
    - threshold : `Int`
        Specifies the point at which a vector has deviated too 
        far from the specified statistical mean or median.
    - median_bool: `Bool`
        If true, specifies that the median should be the turning
        point for the data to be filtered out on. Defaults
        to true. If specified as false, the mean value will be 
        used instead.
        method.
    - m : `Int`
        Defines the number of vectors contributing to the median 
        or mean value of each vector. Defaults to 3, though the
        original implementation mentions that 5 is a good number
        too. Also known as "kernelsize"
    - mask : `Matrix`
        Use to mask out areas of the given matrices to improve
        computation times. Default to an empty matrix.
    Returns:
    --------
    - hu, hv : `Matrices`
            Successfully filtered matrices. New versions of u
            and v.
"""
function localfilt(x, y, u::Matrix{Float32}, v::Matrix{Float32}, 
                        threshold, median_bool=true, m=3, mask=[])
    IN = zeros(eltype(u), size(u))

    dim1 = round(Int32, size(u, 1) + 2 * floor(m/2))
    dim2 = round(Int32, size(u, 2) + 2 * floor(m/2))
    nu = zeros(eltype(u), (dim1, dim2)) * NaN
    nv = zeros(eltype(u), (dim1, dim2)) * NaN
    
    # Transfer over data
    from_cols = round(Int32, floor(m/2) + 1)
    minus_rows = round(Int32, floor(m/2))
    nu[from_cols:end-minus_rows, from_cols:end-minus_rows] = u
    nv[from_cols:end-minus_rows, from_cols:end-minus_rows] = v
    
    INx = zeros(eltype(nu), size(nu))
    INx[from_cols: end - minus_rows, from_cols: end - minus_rows] = IN
    
    U2::Matrix{ComplexF32} = nu .+ im .* nv

    ma, na = size(U2)
    histostd = zeros(ComplexF32, size(nu)) 
    histo = zeros(ComplexF32, size(nu))

    for ii in m - 1:1:na - m + 2
        for jj in m - 1:1:ma - m + 2

            if INx[jj, ii] != 1
                m_floor_two = floor(Int32, m / 2)
                tmp = U2[jj - m_floor_two: jj + m_floor_two,
                        ii - m_floor_two: ii + m_floor_two] 
                tmp[ceil(Int32, m / 2), ceil(Int32, m / 2)] = NaN;

                # Run the appropriate stat depending on method arg.
                usum = median_bool ? im_median_magnitude(tmp[:]) : mean(tmp[:])
                histostd[jj, ii] = im_std(tmp[:])

            else
                usum = NaN; tmp = NaN; histostd[jj, ii] = NaN
            end
            histo[jj, ii] = usum
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

    m_ceil_two = ceil(Int32, m/2)
    m_floor_two = floor(Int32, m/2)
    hu::Matrix{Float32} = nu[m_ceil_two:end - m_floor_two, m_ceil_two:end - m_floor_two]
    hv::Matrix{Float32} = nv[m_ceil_two:end - m_floor_two, m_ceil_two:end - m_floor_two]

    return hu, hv
end

"""
### intpeak
    Interpolates correlation peaks in PIV.
    Parameters:
        x1, y1 : Maximual values in respective directions
        N : Interrogation window size
        Rxm1, Rxp1 : X-max values in matrix R "minus" or "plus" 1.
        Rym1, Ryp1 : Same as above but Y vals.
        R: Matrix resulting from xcorrf2.
    Interpolation uses Gaussian method.

    Original Author:
    Time stamp: 12:32, Apr. 14, 2004.
    Copyright 1998-2004, J. Kristian Sveen, 
    jks@math.uio.no/jks36@damtp.cam.ac.uk
    Dept of Mathmatics, University of Oslo/ 
    DAMTP, Univ. of Cambridge, UK
    Distributed under the GNU general public license.
"""
function intpeak(x1, y1, R, Rxm1, Rxp1, Rym1, Ryp1, N)
    if length(N) == 2
        M = N[1]; N = N[2]
    else
        M = N
    end
    x01 = x1 + ((log(Complex(Rxm1)) - log(Complex(Rxp1))) / ((2 * log(Complex(Rxm1))) - (4 * log(R)) + (2 * log(Complex(Rxp1)))))
    y01 = y1 + ((log(Complex(Rym1)) - log(Complex(Ryp1))) / ((2 * log(Complex(Rym1))) - (4 * log(R)) + (2 * log(Complex(Ryp1)))))
    x0 = x01 - M
    y0 = y01 - N

    x0 = real(x0)
    y0 = real(y0)

    return x0, y0
end


# STATISTICS
"""
## im_median_magnitude
    Take the median of a collection of complex numbers using the absolute magnitude.
    This great function was created by the Julia Community, specifically:
    @PeterSimmon & @mbauman 07/23/2024
    https://discourse.julialang.org/t/median-of-complex-numbers-different-from-matlab-output/117352/5
"""
function im_median_magnitude(collection::AbstractArray{Complex{T}}) where {T}
    i = filter(x -> !isnan(x), collection)
    isempty(i) && return NaN
    n = length(i)
    v = partialsort!(i, div(n+1, 2, RoundDown):div(n+1, 2, RoundUp); by=x -> (abs2(x), angle(x)))
    return sum(v)/length(v)
end

"""
### im_mean
    Find the mean of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_mean(collection)
    if length(collection) < 1
        return NaN
    end
    real_part = mean(real.(collection))
    im_part = mean(imag.(collection))
    return real_part + im_part * im
end

"""
### im_std
    Find the std dev of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
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

"""
        nan_std(collection)

    Compute the standard deviation of a collection, excluding any NaN values.

    Arguments
    ---------
        `collection`: The collection of values.

    Returns
    ---------
        The standard deviation of the collection, excluding NaN values. If the 
        collection is empty or contains only NaN values, NaN is returned.
"""
function nan_std(collection)
    i = filter(x -> !isnan(x), collection)
    if length(i) > 0
        return std(i, corrected=false)
    end
    return NaN
end

"""
        nan_mean(collection)

    Compute the mean of a collection, excluding any NaN values.

    Arguments
    ---------
    - `collection`: A collection of values.

    Returns
    --------
    - The mean of the collection, excluding NaN values. If the collection is 
    empty or contains only NaN values, NaN is returned.
"""
function nan_mean(collection)
    i = filter(x -> !isnan(x), collection)
    if length(i) > 0
        return mean(i)
    end
    return NaN
end


# MAIN
"""
### Main Entry
    Minimal PIV calculation for testing and development.\n
    Run two images through the PIV process, eventually ending in the expected PIV \
    plots.\n
    This should be the only public facing function in this library.\n
    \n**:params:**\n
    A: Matrix containing image data of first frame.\n
    B: Matrix containing image data of second frame.\n
    \n**:returns:**\n
    PIV plots.
"""
function main(A::Matrix{T}, B::Matrix{T}) where {T}
    # Convert the images to matrices of floats
    A = convert(Matrix{Float32}, A)
    B = convert(Matrix{Float32}, B)

    pivwin = 16
    log2pivwin::Int32 = log2(pivwin)
    if log2pivwin - round(log2pivwin) != 0
        error("pivwin must be factor of 2")
    end

    pass_sizes::Vector{Int32} = 2 .^ collect(Int32, 6:-1:log2pivwin)
    push!(pass_sizes, pass_sizes[end]) # Duplicate final element

    # other input params for piv
    dt::Int32 = 1; overlap::Float32 = 0.5; validvec::Int32 = 3
    x, y, u, v, SnR, Pkh = multipassx(A, B, pass_sizes, dt, overlap, validvec)

    # Reject data with too-low signal-to-noise level
    snrthresh::Float32 = 1.3  # OG
    ibad = findall(x -> x < snrthresh, SnR)
    u[ibad] .= NaN
    v[ibad] .= NaN

    # Reject data with too-low correlation peak height
    pkhthresh::Float32 = 0.3
    ibad = findall(x -> x < pkhthresh, Pkh)
    u[ibad] .= NaN
    v[ibad] .= NaN

    # Reject data that disagree strongly with their neighbors in a local window
    u, v = globfilt(u, v)
    
    u_map = heatmap(u, 
                    title = "u [pixels/frame]", 
                    aspect_ratio = :equal, 
                    limits=(0, 200), 
                    xlimits=(0, 385))

    v_map = heatmap(v, 
                    title = "v [pixels/frame]", 
                    aspect_ratio = :equal, 
                    ylimits=(0, 200), 
                    xlimits=(0, 385))
    
    # Display side-by-side
    display(plot(u_map, v_map, layout = (2, 1)))
end


# ------ MAIN SCRIPT ------
function run_piv()

    IM1::String = "juliaPIV/data/im1.jpg"
    IM2::String = "juliaPIV/data/im2.jpg"
    IM1_M::Matrix{Gray{N0f8}} = load(IM1)
    IM2_M::Matrix{Gray{N0f8}} = load(IM2)

    main(IM1_M, IM2_M)
end

timed_run() = @time run_piv()
timed_run()
