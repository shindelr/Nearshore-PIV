"""  
**c = xcorrf2(a,b, pad)**\n
Two-dimensional cross-correlation using Fourier transforms.
    XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
    XCORRF2(A) is the autocorrelation function.
    This routine is functionally equivalent to xcorr2 but usually faster.
    See also XCORR2.
    \n**:params:**
        a: matrix to be compared.\n
        b: matrix to be compared.\n
        pad: Transform and trim result?\n
    \n**:return:** 
        c: A matrix whose values reflect the 2D correlation between every cell
           in a & b.

    Author(s): R. Johnson
    Revision: 1.0   Date: 1995/11/27
"""
function xcorrf2(a, b, pad)
    
end