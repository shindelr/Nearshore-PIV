clear;
function c = xcorrf2(a,b,pad)
% Refactored pad to a default true arg in julia, causing
% the new version to not have control statements accounting
% for whether or not pad is necessary. Could be a mistake. %
if nargin==2
  pad='yes';
end

[ma,na] = size(a);
% disp([ma, na]);
[mb,nb] = size(b);

% make reverse conjugate of one array
b = conj(b(mb:-1:1,nb:-1:1));
% disp(b);
if strcmp(pad,'yes')   % use power of 2 transform lengths
  mf = 2^nextpow2(ma+mb);
  nf = 2^nextpow2(na+nb);
  at = fft2(b,mf,nf);
  bt = fft2(a,mf,nf);
elseif strcmp(pad,'no')
  at = fft2(b);
  bt = fft2(a);
else
  error('Wrong input to XCORRF2');
end

% if pad == 'no' I get an error here with array size
% mismatch?
% multiply transforms then inverse transform
c = ifft2(at.*bt);

% make real output for real input
if ~any(any(imag(a))) && ~any(any(imag(b)))
  c = real(c);
end

% ------------------------------------------------------------
% ~       logical NOT
% any()   returns 1 (true) or 0 (false) looking for the args
% imag()  returns imaginary part of arg
% &&      logical AND
% real()  returns real part of complex arg
% ------------------------------------------------------------

% disp([ma, mb])
% disp([na, nb]);
% disp(mf);
% disp(nf);
% disp(at);
% disp(bt);

if strcmp(pad,'yes');    % trim to standard size
  c(ma+mb:mf,:) = [];
  c(:,na+nb:nf) = [];
elseif strcmp(pad,'no');
  c=(c(1:end-1,1:end-1));
end
end


% ------ TEST ZONE ------
a = [1 2 3 4 5;
  16 2 3 13 2;
  5 11 10 8 9;
  9 7 6 12 4;
  4 14 15 1 7
  ];

b = [1 2 3 4 5;
  6 7 8 9 10;
  11 12 13 14 15;
  16 17 18 19 20;
  21 22 23 24 25
  ];

f = @() xcorrf2(a, b, 'no');
disp(timeit(f));
% disp(xcorrf2(a,b));

% ------ TEST ZONE ------
