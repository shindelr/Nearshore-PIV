function y = m_nanmedian(x)
% MNANMEDIAN NaN protected median value.
%   MNANMEDIAN(X) returns the median treating NaNs as missing values.
%   For vectors, MNANMEDIAN(X) is the median value of the non-NaN
%   elements in X.  For matrices, MNANMEDIAN(X) is a row vector
%   containing the median value of each column, ignoring NaNs.
%
%   See also NANMEAN, NANSTD, NANMIN, NANMAX, NANSUM.

%   by John Peter Acklam, jacklam@math.uio.no
%   

i = ~isnan(x);
if any(i)
  y = median(x(i));
else
  y = NaN;
end

end

og_matrix = [
    0.5488135 + 0.71518937i,  0.60276338 + 0.54488318i,  0.42365480 + 0.64589411i;
    0.43758721 + 0.89177300i,  0.96366276 + 0.38344152i,  0.79172504 + 0.52889492i;
    0.56804456 + 0.92559664i, 0.07103606 + 0.08712930i, 0.02021840 + 0.83261985i
];

nan_matrix = [
    0.5488135 + 0.71518937i,  NaN, 0.42365480 + 0.64589411i;
    0.43758721 + 0.89177300i, 0.96366276 + 0.38344152i, NaN;
    0.56804456 + 0.92559664i, NaN, 0.02021840 + 0.83261985i
];

m_nanmedian(nan_matrix(:))
im_part = [0.71518937, 0.891773, 0.92559664, 0.38344152, 0.64589411, 0.83261985];
real_part = [0.5488, 0.4376, 0.5680, 0.9637, 0.4237, 0.0202];
% disp(median(real_part));