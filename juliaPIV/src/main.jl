using FFTW

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
    \n**:returns:**\n
    x: \n
    y: \n
    datax: \n
    datay: \n
"""
function firstpass(A, B, N, overlap, idx, idy)
    M = N[1]; N = N[2]
    sy, sx = size(A)

    xx_dim1 = ceil(Int64, ((size(A,1)-N) / ((1-overlap) * N))) + 1
    xx_dim2 = ceil(Int64, ((size(A,2)-M) / ((1-overlap) * M))) + 1
    xx = zeros(eltype(A), (xx_dim1, xx_dim2))

    # Some pretty strange code here, but I didn't want to change it just in case
    yy = xx
    datax = xx; datay = xx; IN = zeros(eltype(A), size(A))
    println(((1-overlap) * N), " to ", sy - N + 1)
    cj = 1
    for jj in 1:((1-overlap) * N):sy - N + 1
        # Left off here, range function looks good, weird values though.
        # Same as matlab so no worries.
    end


    # Dummy values
    x=0; y=0;
    return x, y, datax, datay
end




# UTILITIES
"""
### xcorrf2
    Two-dimensional cross-correlation using Fourier transforms.
    XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
    \n**:params:**\n
    A: matrix (2D array) to be compared.\n
    B: matrix ((2D array)) to be compared.\n
    pad: Transform and trim result to optimize the speed of the FFTs. This was faster\
    in matlab. Julia does not handling padding automatically and so required two \
    O(2n) actions to be taken. Default val is false. 
    \n**:return:**\n
    c: A matrix whose values reflect the 2D correlation between every cell \
    in A & B. Matrix is 2D float 64.

    Originally written in Matlab by,\n
    Author(s): R. Johnson\n
    Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(A, B, pad=false)
    # Unpack size() return tuple into appropriate variables
    ma, na = size(A)
    mb, nb = size(B)
    
    # Reverse conjugate
    B = conj(B[mb:-1:1, nb:-1:1])

    # This is A time hog, so pad is default to false.
    if pad
        mf = nextpow(2, ma + mb)  
        nf = nextpow(2, na + nb)

        # Initialize new zero matrix of relevant padding size
        pad_matrix_a = zeros(eltype(A), (mf, nf))
        pad_matrix_b = zeros(eltype(B), (mf, nf))

        # Pad both matrices up to the efficient padding size
        pad_matrix_a[1:size(A,1), 1:size(A,2)] = A[1:size(A,1), 1:size(A,2)]
        pad_matrix_b[1:size(B,1), 1:size(B,2)] = B[1:size(B,1), 1:size(B,2)]

        # Run fft's
        at = fft(pad_matrix_b)
        bt = fft(pad_matrix_a)
    else
        bt = fft(A)
        at = fft(B)
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
    println("===============\nBegin execution\n===============")

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
    println("===============\nExiting now\n===============")
end



# ------ TEST ZONE ------

# A = [1 2 3 4 5;
# 16 2 3 13 2;
# 5 11 10 8 9;  
# 9 7 6 12 4;
# 4 14 15 1 7
# ]

A = [
    38 28 14 42 7 20 38 18 22 10;
    10 23 35 39 23 2 21 1 23 43;
    29 37 1 20 32 11 21 43 24 48;
    26 41 27 15 14 46 50 43 2 36;
    50 6 20 8 38 17 3 24 13 49;
    8 25 1 19 27 46 6 43 7 46;
    34 13 16 35 49 39 3 1 5 41;
    3 28 17 25 43 33 9 35 13 30;
    47 14  7 13 22 39 20 15 44 17;
    46 23 25 24 44 40 28 14 44  0
]

# B = [3 4 5 6 7;
# 8 9 10 11 12;
# 13 14 15 16 17;
# 18 19 20 21 22;
# 23 24 25 26 27
# ]

B = [
    46 34 22 49 7 27 45 28 24 10;
    17 25 37 39 33 6 30 7 32 51;
    35 45  8 21 32 17 27 50 28 50;
    33 46 37 17 14 48 54 45 2 40;
    59 12 26 18 46 26 12 26 19 49;
    11 28  5 25 33 56 9 49 17 48;
    39 14 25 43 53 44 6 11 14 47;
    11 34 17 25 51 43 17 38 21 32;
    53 19 14 23 30 43 20 17 53 24;
    56 28 32 32 47 40 28 23 47 6
]

main(A, B)
# ------ TEST ZONE ------