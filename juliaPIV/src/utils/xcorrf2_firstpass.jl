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
        pad: Transform and trim result? Default val is true. Bool flipped to \
         false if padding is not desired.\n
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




end

# ------ TEST ZONE ------
a = fill(1.0, (5,5))
b = fill(1.0, (3,2))
xcorrf2(a, b)




# ------ TEST ZONE ------