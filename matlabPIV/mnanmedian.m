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

% Robin's code: 
end

test_m = [
  1.0 + 0.0i	1.0 + 0.0i	1.0 + 0.0i;
  1.0 + 0.0i	NaN + 0.0i	0.0 + 0.0i;
  0.0 + 0.0i	0.0 + 0.0i	0.0 + 0.0i
];
m_nanmedian(test_m(:))
