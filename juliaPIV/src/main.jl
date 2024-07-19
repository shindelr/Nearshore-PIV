# Third party modules
# Allows developers to 
using Statistics
using FFTW            # Fast Fourier Transforms library built on C
using Images          # Basic image processing library
using FileIO          # I/O library
using DelimitedFiles  # Write matrices to CSV
using ProgressBars
using Skipper         # Special skipping library to skip NaNs and other things
using Interpolations
using ScatteredInterpolation
using Plots


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
function multipassx(A, B, wins, Dt, overlap, sensit)
    # Convert the images to matrices to floats
    A = convert(Matrix{Float64}, A)
    B = convert(Matrix{Float64}, B)

    sy, sx = size(A)
    total_passes = size(wins, 1)

    # Initial passes are for removing large-scale displacements.  Initialize
    # displacements (datax,datay) to zero
    data_dim_1 = floor(Int64, (sy/(wins[1,1] * (1-overlap))))
    data_dim_2 = floor(Int64, (sx/(wins[1,2] * (1-overlap))))
    datax = zeros(eltype(A), (data_dim_1, data_dim_2))
    datay = copy(datax)
    # for i in 1:total_passes - 1
        # println("Pass ", i, " of ", total_passes )
    i = 1
    x, y, datax, datay = firstpass(A, B, wins[i, :], overlap, datax, datay)
    # TESTING 07/17: Success! First iteration is a perfect match!

    datax, datay = localfilt(x, y, datax, datay, sensit)
    # TESTING 07/18: Down to 10 differences!
    writedlm("tests/juliaOut/multipass_loop/localfilt_datax.csv", datax, ',')

    #     # Not currently working on second iteration?
    #     # datax = naninterp(datax, i)
    #     # datay = naninterp(datay, i)
    #     # TESTING: 60 differences with Multiquadratic() -- Interesting blobs in top zone
    #     # TESTING: 53 differences with InverseMultiquadratic() -- Might be best balance between top zone and center
    #     # TESTING: 53 differences with Gaussian() -- BAD
    #     # TESTING: 43 difference with linear_naninterp() -- Closest to MATLAB
    #     # TESTING: 50 differences with InverseQuadratic() -- Much more uniform extrapolation into the top but more spotty in the center
    #     # TESTING: 52 differences with Polyharmonic() -- weird blob in top left

    #     # OG MATLAB IMPLEMENTATION
        # datax, datay = linear_naninterp(datax, datay)

    #     # TESTING: Down to 11 differences!!

        # datax = floor.(Int, datax)
        # datay = floor.(Int, datay)
        # writedlm("tests/juliaOut/1stpass_linnaninterp_datax.csv", datax, ',')

        # if i != total_passes - 1
            # X, Y, XI, YI = build_grids(wins, overlap, sx, sy, i)

            # STILL NEEDS TESTING, BUT IS NOT CRASHING AT LEAST
            # datax = regular_interp(datax, X, Y, XI, YI)
            # datay = regular_interp(datay, X, Y, XI, YI)
        # end
    # end

    # writedlm("tests/juliaOut/penultimate_datax.csv", datax, ',')
    # writedlm("tests/juliaOut/multipass_loop/penultimate_datay.csv", datay, ',')

    # println("Final Pass")
    # datax = readdlm("tests/juliaOut/multipass_loop/penultimate_datax.csv", ',', Float64)
    # datay = readdlm("tests/juliaOut/multipass_loop/penultimate_datay.csv", ',', Float64)

    # x, y, u, v, SnR, Pkh = finalpass(A, B, wins[end, :], overlap, datax, datay, Dt)


    # Dummy values
    x=0; y=0; u=0; v=0; SnR=0; Pkh=0;
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
    pad: Bool to handle padding of the fast fourier transforms. Default to true,\
    takes extra time but dramatically increases final resolution of the plots.
    \n**:returns:**\n
    x: \n
    y: \n
    datax: \n
    datay: \n
"""
function firstpass(A, B, N, overlap, idx, idy, pad=true)
    M = floor(Int, N[1]); N = floor(Int, N[2])
    if pad
        pad_matrix = pad_for_xcorr( A[1:M, 1:N])
        P = plan_fft(pad_matrix; flags=FFTW.MEASURE)
    else
        P = plan_fft(A[1:M, 1:N])
    end

    # Initializing matrices
    sy, sx = size(A)
    xx_dim1 = ceil(Int, ((size(A,1)-N) / ((1-overlap) * N))) + 1
    xx_dim2 = ceil(Int, ((size(A,2)-M) / ((1-overlap) * M))) + 1
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))
    yy = copy(xx)
    datax = copy(xx)
    datay = copy(xx)
    IN = zeros(Int64, size(A))

    cj = 1
    for jj in 1:((1-overlap) * N):(sy - N + 1)
        ci = 1
        for ii in 1:((1-overlap) * M):(sx - M + 1)
            # Floor correct?
            IN_i_1 = floor(Int64, (jj + N/2))
            IN_i_2 = floor(Int64, (ii + M/2))

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

                # Floor correct?
                C = A[floor(Int, jj):floor(Int, jj+N-1), 
                      floor(Int, ii):floor(Int, ii+M-1)]
                D = B[floor(Int, jj+idy[cj, ci]):floor(Int, jj+N-1+idy[cj, ci]), 
                      floor(Int, ii+idx[cj, ci]):floor(Int, ii+M-1+idx[cj, ci])]
                
                C = C.-mean(C); D = D.-mean(D)
                stad1 = std(C); stad2 = std(D)

                if stad1 == 0
                    stad1 = NaN
                end
                if stad2 == 0
                    stad2 = NaN
                end

                # Call xcorrf2, passing in the FFT plan and normalize result
                if pad
                    R = xcorrf2(C, D, P, true) ./ ( N * M * stad1 * stad2)
                else
                    R = xcorrf2(C, D, P, false) ./ ( N * M * stad1 * stad2)
                end

                # Find position of maximal value of R
                if size(R, 1) == (N - 1)  # I think checking for second to last pass
                    max = maximum(R)
                    max_coords = findall(x -> x == max, R)
                else
                    subset = R[Int(floor(.5*N+2)):Int(floor(1.5*N-3)), 
                                Int(floor(.5*M+2)):Int(floor(1.5*M-3))]
                    # max = maximum(R)
                    max = maximum(subset)
                    max_coords = findall(x -> x == max, subset)
                    max_coords = [(i[1] + Int(0.5*N+1), i[2] + Int(0.5*M+1)) for i in max_coords]
                end

                # Handle a vector that has multiple maximum coordinates.
                # Sum the product of each x and y indice with its own indice within
                # the max_coords vector.
                if length(max_coords) > 1

                    max_x1 = round(
                        Int, sum([c[2] * i for (c, i) in enumerate(max_coords)]) / 
                        sum([c[2] for c in max_coords]))
                    max_y1 = round(
                        Int, sum([c[1] * i for (c, i) in enumerate(max_coords)]) / 
                        sum([c[1] for c in max_coords]))

                # Handle empty max_coords vector.
                elseif isempty(max_coords)
                    idx[cj, ci] = NaN
                    idy[cj, ci] = NaN
                    max_x1 = NaN
                    max_y1 = NaN
                
                # Otherwise, unpack into max coordinates to be used in code below.
                else
                    max_y1, max_x1 = max_coords[1]
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

# HALFWAY FINISHED
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
function finalpass(A, B, N, ol, idx, idy, Dt, pad=true)
    # Set up
    if length(N) == 1
        M = N
    else
        M = N[1]; N = N[2]
    end

    if pad
        pad_matrix = pad_for_xcorr(A[1:convert(Int, M), 1:convert(Int, N)])
        P = plan_fft(pad_matrix; flags=FFTW.MEASURE)
    else
        P = plan_fft(A[1:M, 1:N])
    end

    cj = 1
    sy, sx = size(A)
    dim_1 = ceil(Int, (sy - N) / ((1 - ol) * N)) + 1
    dim_2 = ceil(Int, (sx - M) / ((1 - ol) * M)) + 1

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
                floor(Int, jj + idy[cj, ci]):floor(Int, jj + N - 1 + idy[cj, ci]),
                floor(Int,ii + idx[cj, ci]): floor(Int, ii + M - 1 + idx[cj, ci])
                ]
            E = A[floor(Int, jj):floor(Int, jj + N - 1),
                  floor(Int, ii):floor(Int, ii + M - 1)]

            stad1 = std(E); stad2 = std(D2)

            if stad1 == 0
                stad1 = 1
            end
            if stad2 == 0
                stad2 = 1
            end

            E = E.-mean(E); F = D2.-mean(D2)

            # Cross correlate and FFT
            if pad
                R = xcorrf2(E, F, P, true) ./ ( N * M * stad1 * stad2)
            else
                R = xcorrf2(E, F, P, false) ./ ( N * M * stad1 * stad2)
            end

            if !any(isnan.(R)) & !all(x -> x == 0, R)

                # Find position of maximal value of R
                if size(R, 1) == (N - 1)
                    max = maximum(R)
                    max_coords = findall(x -> x == max, R)
                else
                    subset = R[Int(floor(.5 * N + 2)):Int(floor(1.5 * N - 3)), 
                                Int(floor(.5 * M + 2)):Int(floor(1.5 * M - 3))]
                    # max = maximum(R)
                    max = maximum(subset)
                    max_coords = findall(x -> x == max, subset)
                    max_coords = [(i[1] + Int(0.5*N+1), i[2] + Int(0.5*M+1)) for i in max_coords]
                end

                if length(max_coords) == 0
                    @show max max_coords
                end


                # Handle a vector that has multiple maximum coordinates.
                # Sum the product of each x and y indice with its own indice within
                # the max_coords vector.
                if length(max_coords) > 1
                    max_x1 = round(Int, 
                            sum([c[2]^2 for (c, i) in enumerate(max_coords)]) ./
                            sum([c[2] for c in max_coords]))
                    max_y1 = round(Int, 
                            sum([c[1]^2 for (c, i) in enumerate(max_coords)]) ./
                            sum([c[1] for c in max_coords]))
                end
                
                # Unpack cartesian index type. Only a handful iterations here
                if length(max_coords) == 1
                    for (i, j) in max_coords
                        max_y1 = i; max_x1 = j
                    end
                end
                
                # Some kind of manual adjustment?
                # if max_x1 == 1
                #     max_x1 = 2
                # end
                # if max_y1 == 1
                #     max_y1 = 2
                # end




            else
            end


            ci += 1
        end
        cj += 1
    end
    
    println("Leaving Final Pass")
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
function pad_for_xcorr(trunc_matrix)
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
    pad: Transform and trim result to optimize the speed of the FFTs and\
    increase the resolution of the resulting PIV plots. This was faster in \
    matlab. Default val is true and the preoptimized FFT plan is used to offset\
    the time issues. \n
    plan: Callable function representing a pre-planned, omptimized FFT. Created \
    using the dimensions of matrices A & B.
    \n**:return:**\n
    c: A matrix whose values reflect the 2D correlation between every cell \
    in A & B. Matrix is 2D float 64.

    Originally written in Matlab by,\n
    Author(s): R. Johnson\n
    Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(A, B, plan, pad=true)
    # Unpack size() return tuple into appropriate variables
    ma, na = size(A)
    mb, nb = size(B)
    
    # Reverse conjugate
    B = conj(B[mb:-1:1, nb:-1:1])

    # This is a time hog, but increases resolution of final plots by quite a bit
    # Room for improvement here. I could pass in the original padded matrix somehow.
    if pad
        pad_matrix_a = pad_for_xcorr(A)
        pad_matrix_b = pad_for_xcorr(B)

        # Transfer data from og matrix to optimized sized ones
        pad_matrix_a[1:size(A,1), 1:size(A,2)] = A[1:size(A,1), 1:size(A,2)]
        pad_matrix_b[1:size(B,1), 1:size(B,2)] = B[1:size(B,1), 1:size(B,2)]

        # Runs optimized FFTs using the plan
        at = plan * pad_matrix_b
        bt = plan * pad_matrix_a
    else
        bt = plan * A
        at = plan * B
    end

    # Mult transforms and invert
    mult_at_bt = at.*bt
    c = ifft(mult_at_bt)

    # Make all real
    c = real(c)

    # Trim
    if pad
        rows = ma + mb
        cols = na + nb
        # Keep everything from the first index to the index where padding began
        c = c[1:rows - 1, 1:cols - 1]
    else
        c = c[1:end-1, 1:end-1]
    end
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
function linear_naninterp(u, v)
    coords = findall(x->isnan(x), u)
    numm = length(coords)
    dy,dx = size(u)
    lp = 1; tel = 1

    # counter = 0
    # Now sort the NaN's after how many neighbors they have that are
    # physical values. Then we first interpolate those that have 8
    # neighbors, followed by 7, 6, 5, 4, 3, 2 and 1

    while !isempty(coords)
        nei = zeros(Int64, length(coords), 3)

        # writedlm("tests/juliaOut/Jtest_py_px.csv", 
        #         [(coords[i][1],coords[i][2]) for i in eachindex(coords)], 
        #         ',')

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

        # !!! TESTING: Strange results, looks the same by hand but off by 24
        # rows. Seems like the nan count is off by one in some cases?
        # writedlm("tests/juliaOut/Jtest_presort_NEI.csv", nei, ',')
        
        # Sort NEI by row to interpolate vectors with fewest spurious neighbors.
        nei = sortslices(nei, dims=1, lt=Base.isgreater)
        # writedlm("tests/juliaOut/Jtest_postsort_NEI.csv", nei, ',')
  
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

            temp_u = u[py[j] - 1 + cory1:py[j] + 1 + cory2, 
                    px[j] - 1 + corx1:px[j] + 1 + corx2] 
            temp_v = v[py[j] - 1 + cory1:py[j] + 1 + cory2, 
                    px[j] - 1 + corx1:px[j] + 1 + corx2] 
            
            mean_prep_u = collect(Skipper.skip(x -> isnan(x), temp_u[:]))
            mean_prep_v = collect(Skipper.skip(x -> isnan(x), temp_v[:]))
            u[py[j], px[j]] = mean(mean_prep_u)
            v[py[j], px[j]] = mean(mean_prep_v)
            
            if lp > numm
                u[py[j], px[j]] = 0
                v[py[j], px[j]] = 0
            end
        end
        tt = length(py)
        coords = findall(x->isnan(x), u)  # Might not be necessary
        lp += 1
        # counter += 1
        # println(counter)
        # display(coords)
    end
    return u, v
end

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

function regular_interp(samples, xs, ys, XI, YI)
    itp = Interpolations.interpolate((ys, xs), samples, Gridded(Linear()))
    extp = Interpolations.extrapolate(itp, Line())

    itp_results = zeros(Float64, (length(YI), length(XI)))

    for (mi, yi) in enumerate(YI)
        for (ni, xi) in enumerate(XI)
            # Interpolate the interior of the matrix
            if (33 < xi < 3041) && (33 < yi < 2017)
                itp_results[mi, ni] = itp(yi, xi)
            else
            # Extrapolate the exterior
                itp_results[mi, ni] = extp(yi, xi)
            end
        end
    end

    # return round.(Int, itp_results)
    # Returning float to stay compatible with the rest of the program?
    return itp_results

end

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
function localfilt(x, y, u, v, threshold, median_bool=true, m=3, mask=[])
    method =  median_bool ? "median" : "mean"
    IN = zeros(eltype(u), size(u))

    dim1 = round(Int, size(u, 1) + 2 * floor(m/2))
    dim2 = round(Int, size(u, 2) + 2 * floor(m/2))
    nu = zeros(eltype(u), (dim1, dim2)) * NaN
    nv = zeros(eltype(u), (dim1, dim2)) * NaN
    
    # Transfer over data
    from_cols = round(Int, floor(m/2) + 1)
    minus_rows = round(Int, floor(m/2))
    nu[from_cols:end-minus_rows, from_cols:end-minus_rows] = u
    nv[from_cols:end-minus_rows, from_cols:end-minus_rows] = v
    # TESTING 07/28: Success! NV/NU are both equivalent to matlab
    
    INx = zeros(eltype(nu), size(nu))
    INx[from_cols: end - minus_rows, from_cols: end - minus_rows] = IN
    # Testing 07/28: Success! INx equivalent to matlab
    
    # Not sure any of these vars are used.
    # prev = isnan.(nu)
    # previndex = findall(prev)
    # teller = true

    U2 = nu .+ im .* nv
    # Testing 07/28: Sucess! Matlab equivalency. 
    # writedlm("tests/juliaOut/first_localfilt/U2.csv", U2, ',')

    ma, na = size(U2)
    histostd = zeros(ComplexF64, size(nu)) 
    histo = zeros(ComplexF64, size(nu))

    # Not sure these two are used either
    # hista = zeros(eltype(nu), size(nu)) 
    # histastd = zeros(eltype(nu), size(nu)) 

    iter = ProgressBar(m - 1:na - m + 2)
    for _ in iter  # Looks gnar, but just a bar!
        for ii in m - 1:1:na - m + 2
            for jj in m - 1:1:ma - m + 2

                if INx[jj, ii] != 1
                    m_floor_two = floor(Int, m / 2)
                    tmp = U2[jj - m_floor_two: jj + m_floor_two,
                            ii - m_floor_two: ii + m_floor_two] 
                    tmp[ceil(Int, m / 2), ceil(Int, m / 2)] = NaN;

                    # Create a collection of all elements without NaN values
                    usum_prep = collect(Skipper.skip(x -> isnan(x), tmp[:]))

                    # Run the appropriate stat depending on method arg.
                    usum = median_bool ? im_median(tmp[:]) : mean(tmp[:])
                    histostd[jj, ii] = im_std(usum_prep)

                else
                    usum = NaN; tmp = NaN; histostd[jj, ii] = NaN
                end
                histo[jj, ii] = usum
            end
        end 
        set_description(iter, "Local $method filter running: ")
    end

    # TESTING 07/08: histostd matrices are equivalent to 11 decimal points, then 
    #          matlab rounds off and Julia continues on for a few more digits
    writedlm("tests/juliaOut/first_localfilt/histostd.csv", histostd, ',')

    # TESTING 07/08: histo matrices are off by 50 rows! I think this could be
    #           where we're getting messed up. Needs more investigation.
    writedlm("tests/juliaOut/first_localfilt/histo.csv", histo, ',')
    
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

    # TESTING: Showing 13 differences, nu & nv are directly related to datax
    #  and datay. So maybe this problem stems from firstpass()

    writedlm("tests/juliaOut/first_localfilt/nu.csv", nu, ',')
    writedlm("tests/juliaOut/first_localfilt/nv.csv", nv, ',')

    # Skipped print statement about how many vectors were filtered.
    # Skpped checking for 'interp' arg, because the actual program wasn't using
    # it. We call naninterp explicitly right after this function.

    m_ceil_two = ceil(Int, m/2)
    m_floor_two = floor(Int, m/2)
    hu = nu[m_ceil_two:end - m_floor_two, m_ceil_two:end - m_floor_two]
    hv = nv[m_ceil_two:end - m_floor_two, m_ceil_two:end - m_floor_two]
    return hu, hv
end


# COMPLEX NUMBER STATISTICS
"""
### im_median
    Find the median of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_median(collection)
    i = filter(x -> !isnan(x), collection)

    if length(i) > 0
        real_part = median(real(i))
        im_part = median(imag(i))
        return real_part + im_part * im
    end

    # If all NaN
    return NaN
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
### im_median
    Find the std dev of the argued collection of complex numbers.
    If the collection is empty, returns NaN.
"""
function im_std(collection)
    if length(collection) < 1
        return NaN
    end
    real_part = std(real(collection), corrected=false)
    im_part = std(imag(collection), corrected=false)
    return real_part + im_part * im
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
function main(A, B)
    println("===============\nBegin execution\n===============\n")

    pivwin = 16
    log2pivwin = log2(pivwin)
    if log2pivwin - round(log2pivwin) != 0
        error("pivwin must be factor of 2")
    end

    pass_sizes = 2 .^ collect(6:-1:log2pivwin)
    push!(pass_sizes, pass_sizes[end]) # Duplicate final element
    pass_sizes = [pass_sizes pass_sizes] # Convert vector to matrix

    # other input params for piv
    dt = 1; overlap = 0.5; validvec = 3
    @time x, y, u, v, SnR, Pkh = multipassx(A, B, pass_sizes, dt, overlap, validvec)

    println("\n===============\nExiting now\n===============")
end


# ------ TEST ZONE ------
# im1_path = "../data/im1.jpg"
im1_path = "juliaPIV/data/im1.jpg"
im2_path = "juliaPIV/data/im2.jpg"
A = load(im1_path)
B = load(im2_path)

main(A, B)
# ------ TEST ZONE ------