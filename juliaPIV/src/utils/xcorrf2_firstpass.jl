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
        println("MA: ", ma, "\nNA: ", na)
        println("MB: ", mb, "\nNB: ", nb)

        # Make reverse conjugate of array B
        rev_b = b[mb:-1:1, nb:-1:1] 
        b = conj(rev_b);
        println("Reversed Conjugate B:\n", b)
        # Either transform len with pow2 or skip straight to fourier transforms.
        if pad
            # Slightly different behavior of nextpow() eliminates the need for 2^
            mf = nextpow(2, ma+mb)  
            # println("MF: ", mf)
            nf = nextpow(2, na+nb)
            # println("NF: ", nf)
            at = fft(b, [mf, nf])  # Opportunity for more efficiency using plan_fft()?
            # LEFT OFF HERE. ISSUE WITH FFT()

        end




end

# ------ TEST ZONE ------
a = [16 2 3 13;
     5 11 10 8;  
     9 7 6 12;
     4 14 15 1]
b = [1 2 3; 4 5 6; 7 8 9]
xcorrf2(a, b)

# ------ TEST ZONE ------