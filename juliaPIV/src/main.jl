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

    

    x=0; y=0; u=0; v=0; SnR=0; Pkh=0;
    return x, y, u, v, SnR, Pkh

end


"""
TODO:
"""
function firstpass(A, B, N, overlap, idx, idy)
    

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
    sy,sz = size(A)
    x, y, u, v, SnR, Pkh = multipassx(A, B, pass_sizes, dt, overlap, validvec)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!! Profile stuff goes here !!!!!!!!!!!!!
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    

end



# ------ TEST ZONE ------

A = [1 2 3 4 5;
16 2 3 13 2;
5 11 10 8 9;  
9 7 6 12 4;
4 14 15 1 7
]

B = [3 4 5 6 7;
8 9 10 11 12;
13 14 15 16 17;
18 19 20 21 22;
23 24 25 26 27
]

main(A, B)
# ------ TEST ZONE ------