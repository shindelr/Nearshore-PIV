using FFTW
using Statistics

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

    for i in 1:iter-1
        println("Iter ", i, " of ", iter )
        # PIV proper
        x, y, datax, datay = firstpass(A, B, wins[i, :], overlap, datax, datay)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!! Complete multipassx below here !!!!!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end

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
    # Plan_fft using A, M, N here!
    # Chose to use paitient because it finds the perfect optimiazation for the
    # given matrix. It takes a couple of seconds at run time, but then the FFTs
    # are all run at that optimized speed, supposedly making up for the loss.
    
    if pad
        pad_matrix = pad_for_xcorr( A[1:M, 1:N])
        P = plan_fft(pad_matrix; flags=FFTW.PATIENT)
    else
        P = plan_fft(A[1:M, 1:N])
    end


    sy, sx = size(A)
    xx_dim1 = ceil(Int, ((size(A,1)-N) / ((1-overlap) * N))) + 1
    xx_dim2 = ceil(Int, ((size(A,2)-M) / ((1-overlap) * M))) + 1
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))

    # Some pretty strange code here, but I didn't want to change it just in case
    yy = xx
    datax = xx; datay = xx; 
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
                
                # Have to broadcast "." over each element
                C = C.-mean(C); D = D.-mean(D)
                stad1 = std(C); stad2 = std(D)  # Might need to vec() these

                # To apply weight function, uncomment below:
                # C = C.*W; D = D.*W
                
                if stad1 == 0
                    stad1 = NaN
                end
                if stad2 == 0
                    stad2 == NaN
                end

                # Call xcorrf2, passing in the FFT plan and normalize result
                if pad
                    R = xcorrf2(C, D, P, true) ./ ( N* M * stad1 * stad2)
                else
                    R = xcorrf2(C, D, P, false) ./ ( N* M * stad1 * stad2)
                end

                
                ci=ci+1;  # This placement will need to be adjusted
            end
            cj += 1
        end
    end


    # Dummy values
    x=0; y=0;
    return x, y, datax, datay
end


# UTILITIES
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
    return zeros(eltype(A), (mf, mf)) 
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

# A = [
#     38 28 14 42 7 20 38 18 22 10;
#     10 23 35 39 23 2 21 1 23 43;
#     29 37 1 20 32 11 21 43 24 48;
#     26 41 27 15 14 46 50 43 2 36;
#     50 6 20 8 38 17 3 24 13 49;
#     8 25 1 19 27 46 6 43 7 46;
#     34 13 16 35 49 39 3 1 5 41;
#     3 28 17 25 43 33 9 35 13 30;
#     47 14  7 13 22 39 20 15 44 17;
#     46 23 25 24 44 40 28 14 44  0
# ]

# Generate a matrix of random values
A = rand(2048, 3072)
B = A .+ 2

# B = [
#     46 34 22 49 7 27 45 28 24 10;
#     17 25 37 39 33 6 30 7 32 51;
#     35 45  8 21 32 17 27 50 28 50;
#     33 46 37 17 14 48 54 45 2 40;
#     59 12 26 18 46 26 12 26 19 49;
#     11 28  5 25 33 56 9 49 17 48;
#     39 14 25 43 53 44 6 11 14 47;
#     11 34 17 25 51 43 17 38 21 32;
#     53 19 14 23 30 43 20 17 53 24;
#     56 28 32 32 47 40 28 23 47 6
# ]

@time main(A, B)
# ------ TEST ZONE ------