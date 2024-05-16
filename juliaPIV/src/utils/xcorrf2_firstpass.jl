using FFTW
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
pad: Transform and trim result? Default val is true. Pass false if padding is \
not desired.\n
\n**:return:**\n
c: A matrix whose values reflect the 2D correlation between every cell \
in a & b.

Originally written in Matlab by,
Author(s): R. Johnson
Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(a, b, pad=true)
    # Opportunity for more efficiency using plan_fft()?
    
    # Unpack size() return tuple into appropriate variables
    ma, na = size(a)
    mb, nb = size(b)
    
    # Reverse conjugate
    b = conj(b[mb:-1:1, nb:-1:1])

    # Either transform len with pow2 or skip straight to fourier transforms.
    if pad
        # Slightly different behavior of nextpow() eliminates the need for 2^ maybe?
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

    # Pad = False, don't do it
    else
        bt = fft(a)
        at = fft(b)  # Opportunity for more efficiency using plan_fft()?
    end

    # Mult transforms and invert

    # pad = false generates a mismatched size error.
    mult_at_bt = at.*bt

    # Seems to work, but I get a ton of imaginary junk where matlab gets 00.00's
    # I don't think it matters since we trim it off anyways.
    c = ifft(mult_at_bt)

    # I'm not sure this is necessary? I mean won't it be at least O(n) to check
    # both arrays for any imaginary numbers? Seems like a real time hog, and
    # real() runs happily on both complex and real nums.
    # Real out for real in
    # if !any(imag.(a)) && !any(imag.(b))
    #     c = real(c)
    # end

    # Make all real
    real(c)

    println("ma: ", ma,"\nna: ", na,"\nmb: ", mb, "\nnb: ", nb, "\nmf: ", mf, "\nnf: ",nf)

    # Trim
    if pad
        rows = ma + mb
        cols = na + nb

        # This seems to work for the matrices I've tested it on. Needs more 
        # testing though! 
        # We just keep everything from the first index to the index where we
        # started padding things.
    
        c = c[1:rows - 1, 1:cols - 1]
    else
        c = c[1:end-1, 1:end-1]
    end
end

# ------ TEST ZONE ------
a = [16 2 3 13 2;
     5 11 10 8 9;  
     9 7 6 12 4;
     4 14 15 1 7
     ]

# a = [400.0+0.0im 434.0+0.0im;         
# 491.0+0.0im 795.0+0.0im]
# b = [400.0+0.0im 434.0+0.0im;         
# 491.0+0.0im 795.0+0.0im]

b = [1 2 3 4 5;
    6 7 8 9 10;
    11 12 13 14 15;
    16 17 18 19 20;
    21 22 23 24 25
    ]

@time xcorrf2(a, b)

# ------ TEST ZONE ------