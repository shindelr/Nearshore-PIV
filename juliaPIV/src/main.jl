using FFTW

"""
TODO:
"""
function firstpass(A, B, N, overlap, idx, idy)
    

end






"""
**c = xcorrf2(a,b, pad)**\n
Two-dimensional cross-correlation using Fourier transforms.
XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
XCORRF2(A) is the autocorrelation function.
This routine is functionally equivalent to xcorr2 but usually faster.
See also XCORR2.
\n**:params:**\n
a: matrix (2D array) to be compared.\n
b: matrix ((2D array)) to be compared.\n
pad: Transform and trim result to optimize the speed of the FFTs. This was faster\
 in matlab. Julia does not handling padding automatically and so required two \
 O(2n) actions to be taken. Default val is false. 
\n**:return:**\n
c: A matrix whose values reflect the 2D correlation between every cell \
in a & b. matrix is 2D float 64.

Originally written in Matlab by,\n
Author(s): R. Johnson\n
Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(a, b, pad=false)
    # Unpack size() return tuple into appropriate variables
    ma, na = size(a)
    mb, nb = size(b)
    
    # Reverse conjugate
    b = conj(b[mb:-1:1, nb:-1:1])

    # This is a time hog, so pad is default to false.
    if pad
        mf = nextpow(2, ma + mb)  
        nf = nextpow(2, na + nb)

        # Initialize new zero matrix of relevant padding size
        pad_matrix_a = zeros(eltype(a), (mf, nf))
        pad_matrix_b = zeros(eltype(b), (mf, nf))

        # Pad both matrices up to the efficient padding size
        pad_matrix_a[1:size(a,1), 1:size(a,2)] = a[1:size(a,1), 1:size(a,2)]
        pad_matrix_b[1:size(b,1), 1:size(b,2)] = b[1:size(b,1), 1:size(b,2)]

        # Run fft's
        at = fft(pad_matrix_b)
        bt = fft(pad_matrix_a)
    else
        bt = fft(a)
        at = fft(b)
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
end


"""
Minimal PIV calculation for testing and development.
"""
function main()

    # !!!!!!!!!!Replace these test matrices with the loaded image data!!!!!!!!
    a = [1 2 3 4 5;
        16 2 3 13 2;
        5 11 10 8 9;  
        9 7 6 12 4;
        4 14 15 1 7
        ]

    b = [3 4 5 6 7;
        8 9 10 11 12;
        13 14 15 16 17;
        18 19 20 21 22;
        23 24 25 26 27
        ]
    
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

    # !!!!!!!!!! Profile stuff goes here !!!!!!!!!!!!!

    



end



# ------ TEST ZONE ------
main()
# ------ TEST ZONE ------