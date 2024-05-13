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

Author(s): R. Johnson
Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(a, b, pad=true)
    
    # Unpack size() return tuple into appropriate variables
    ma, na = size(a)
    mb, nb = size(b)
    
    # Make reverse conjugate of array B
    rev_b = b[mb:-1:1, nb:-1:1] 
    b = conj(rev_b);

    # Either transform len with pow2 or skip straight to fourier transforms.
    if pad
        # Slightly different behavior of nextpow() eliminates the need for 2^ maybe?
        mf = nextpow(2, ma+mb)  
        nf = nextpow(2, na+nb)

        # Manual padding? Concat not currently working bc diff sized matrices
        b_pad = [b; zeros(eltype(b), (mf, nf))]
        println("B_pad: ", b_pad)

        at = fft(b)  # Opportunity for more efficiency using plan_fft()?
        bt = fft(a)

        # There may be a big time issue with fft here, can't figure out
        # how to add the padding.
        
    # Pad = False, don't do it
    else
        bt = fft(a)
        at = fft(b)  # Opportunity for more efficiency using plan_fft()?
    end
    
    # Can we do these transforms in place?
end

# ------ TEST ZONE ------
a = [16 2 3 13 2;
     5 11 10 8 9;  
     9 7 6 12 4;
     4 14 15 1 7
     ]

b = [1 2 3 4 5;
    6 7 8 9 10;
    11 12 13 14 15;
    16 17 18 19 20;
    21 22 23 24 25
    ]

# Get dims for padding
mb, nb = size(b)
ma, na = size(a)

println("mb: ", mb, "\nnb: ", nb)
# Reverse conjugate
b = conj(b[mb:-1:1, nb:-1:1])
println("Reverse Conjugate: ")
b  # Print b

# Get padding dimensions, next power of 2 to optimize fft
mf = nextpow(2, ma + mb)
nf = nextpow(2, na + nb)
println("mf: ", mf, "\nnf: ", nf)

# Initialize new zero matrix of relevant padding size
pad_matrix_a = zeros(eltype(a), (mf, nf))
pad_matrix_b = zeros(eltype(b), (mf, nf))

# Pad both matrices up to the efficient padding size
pad_matrix_a[1:size(a,1), 1:size(a,2)] = a[1:size(a,1), 1:size(a,2)]
pad_matrix_b[1:size(b,1), 1:size(b,2)] = b[1:size(b,1), 1:size(b,2)]

# Holy moly do the fft!
at = fft(pad_matrix_b)
bt = fft(pad_matrix_b)

# ------ TEST ZONE ------