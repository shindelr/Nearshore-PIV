# Third party modules
# Allows developers to 
using Statistics
using FFTW            # Fast Fourier Transforms library built on C
using Images          # Basic image processing library
using FileIO          # I/O library
using DelimitedFiles  # Write matrices to CSV
using ProgressBars
using Skipper         # Special skipping library to skip NaNs and other things
# using Interpolations
using ScatteredInterpolation
using Plots

# User defined modules
include("./filtering.jl")
using .PIVFilters: localfilt #linear_naninterp

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
    # Convert A, B to floats
    A = convert(Matrix{Float64}, A)
    B = convert(Matrix{Float64}, B)

    sy, sx = size(A)
    iter = size(wins, 1)

    # Initial passes are for removing large-scale displacements.  Initialize
    # displacements (datax,datay) to zero
    data_dim_1 = floor(Int64, (sy/(wins[1,1] * (1-overlap))))
    data_dim_2 = floor(Int64, (sx/(wins[1,2] * (1-overlap)))) 
    datax = zeros(eltype(A), (data_dim_1, data_dim_2))
    datay = zeros(eltype(A), (data_dim_1, data_dim_2))

    # Disabled loop for testing
    # Getting slight errors on datay and datax. Last test showed 5 differences.
        # for i in 1:iter-1
            # println("Iter ", i, " of ", iter )
            # PIV proper
            i = 1

            # writedlm("tests/juliaOut/jtest_DATAX.csv", datax, ',')

            x, y, datax, datay = firstpass(A, B, wins[i, :], overlap, datax, datay)
            
            datax, datay = localfilt(x, y, datax, datay, sensit)
            # TESTING: 13 differences, same as detected prior.
            
            datax, datay = linear_naninterp(datax, datay)
            # TESTING: 43 differences. We're going to need to debug firstpass
            # to see why datax, datay are experiencing differences.
            
            datax = floor.(Int, datax)
            datay = floor.(Int, datay)
            """TESTING: 12-15 differences almost all in the same columns and rows
            as one another. But after flooring the matrices to integers, the
            43 differences above smoothed out, could've just been differences
            in data representation. This problem must stem back to firstpass."""
            # writedlm("tests/juliaOut/jtest_INTERPDATAX.csv", datax, ',')

            # Different process for the final pass
            if i != iter - 1
                next_win_x = wins[i + 1, 1]
                next_win_y = wins[i + 1, 2]

                # Final window size is duplicated, so check for equality.
                if wins[i, 1] != next_win_x
                    # First time: X = [1:32:3009] .+ 32
                    # X = 33, 65, ..., 3041   size = 95
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

                # writedlm("tests/juliaOut/jtest_X.csv", round.(Int,X'), ',')

                # Creating new NaN matrices to interpolate into
                # m = length(YI)
                # n = length(XI)
                # itp_datax = fill(NaN, m+2, n+2)
                # itp_datay = fill(NaN, m+2, n+2)

                # Build interpolate func, layer onto datax, copy into NaN matrix
                # itp_x = interpolate((Y, X), datax, Gridded(Linear()))
                # datax = [itp_x(yi, xi) for yi in YI, xi in XI]
                
                # itp_y = interpolate((Y, X), datay, Gridded(Linear()))
                # datay = [itp_y(yi, xi) for yi in YI, xi in XI]


                
                # itp_datax[2:end-1, 2:end-1] = round.(Int, datax)
                # itp_datay[2:end-1, 2:end-1] = round.(Int, datay)
                
                # writedlm("tests/juliaOut/jtest_DATAX.csv", datax, ',')

                # Interpolate out the NaN's we put in
                # datax, datay = linear_naninterp(itp_datax, itp_datay)


                
                # datax = round.(Int, datax)
                # datay = round.(Int, datay)

            end
        # end

    # Dummy values
    x=0; y=0; u=0; v=0; SnR=0; Pkh=0;
    return x, y, u, v, SnR, Pkh

end

# Should break up this func into multiple functions. It's huge!
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

    sy, sx = size(A)
    xx_dim1 = ceil(Int, ((size(A,1)-N) / ((1-overlap) * N))) + 1
    xx_dim2 = ceil(Int, ((size(A,2)-M) / ((1-overlap) * M))) + 1
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))
    yy = zeros(eltype(A), (xx_dim1, xx_dim2))
    datax = zeros(eltype(A), (xx_dim1, xx_dim2))
    datay = zeros(eltype(A), (xx_dim1, xx_dim2))
    IN = zeros(Int64, size(A))

    cj = 1
    for jj in 1:((1-overlap) * N):(sy - N + 1)
        ci = 1
        for ii in 1:((1-overlap) * M):(sx - M + 1)
            # Using floor until I have more information! Could be a problem.
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
                elseif (jj + idy[cj, ci]) > (sy - N+1)
                    idy[cj, ci] = sy - N + 1 - jj
                end
                if (ii + idx[cj, ci]) < 1
                    idx[cj, ci] = 1 - ii
                elseif (ii + idx[cj, ci]) > (sx - M + 1)
                    idx[cj, ci] = sx - M + 1 - ii
                end

                # Using floor until I have more information! Could be a problem.
                C = A[floor(Int, jj):floor(Int, jj+N-1), 
                      floor(Int, ii):floor(Int, ii+M-1)]
                D = B[floor(Int, jj+idy[cj, ci]):floor(Int, jj+N-1+idy[cj, ci]), 
                      floor(Int, ii+idx[cj, ci]):floor(Int, ii+M-1+idx[cj, ci])]
                
                C = C.-mean(C); D = D.-mean(D)
                stad1 = std(C); stad2 = std(D)  # Might need to vec() these

                if stad1 == 0
                    stad1 = NaN
                end
                if stad2 == 0
                    stad2 == NaN
                end

                # Call xcorrf2, passing in the FFT plan and normalize result
                if pad
                    R = xcorrf2(C, D, P, true) ./ ( N * M * stad1 * stad2)
                else
                    R = xcorrf2(C, D, P, false) ./ ( N * M * stad1 * stad2)
                end
                

            """
                This was super tricky! Matlab and Julia have quite different 
                implementations in this region. 
                Julia has this findall function with a really nifty arrow syntax, 
                but using it returns a strange vector of CartesianIndex objects that 
                are pretty tough to work with. To handle it, I had to use some weird 
                tuple unpacking operations and list comprehension that matlab didn't
                need. I wonder how true Julia users would have done this section.
            """
                
                # Find position of maximal value of R
                if size(R, 1) == (N - 1)
                    max = maximum(R)
                    max_coords = findall(x -> x == max, R)
                else
                    subset = R[Int(floor(.5*N+2)):Int(floor(1.5*N-3)), 
                                Int(floor(.5*M+2)):Int(floor(1.5*M-3))]
                    max = maximum(R)
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

# FILTERS
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

function naninterp(u, v)
    coords = findall(x->isnan(x), u)


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
    x, y, u, v, SnR, Pkh = multipassx(A, B, pass_sizes, dt, overlap, validvec)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!!Profile stuff goes here and rest of main() lol !!!!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    println("\n===============\nExiting now\n===============")
end


# ------ TEST ZONE ------
# im1_path = "../data/im1.jpg"
im1_path = "juliaPIV/data/im1.jpg"
im2_path = "juliaPIV/data/im2.jpg"
A = load(im1_path)
B = load(im2_path)

@time main(A, B)
# ------ TEST ZONE ------