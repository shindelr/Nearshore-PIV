clear;
function c = xcorrf2(a,b,pad)
% Refactored pad to a default true arg in julia, causing 
% the new version to not have control statements accounting
% for whether or not pad is necessary. Could be a mistake. %
% disp(a)
% disp(b)
if nargin==2  
    pad='yes';
end

  [ma,na] = size(a);
  disp([ma, na]);
  [mb,nb] = size(b);
  disp([mb, nb]);

  % make reverse conjugate of one array
  b = conj(b(mb:-1:1,nb:-1:1));
  disp(b);
  if strcmp(pad,'yes')   % use power of 2 transform lengths
    mf = 2^nextpow2(ma+mb);
    % disp(mf);
    nf = 2^nextpow2(na+nb);
    % disp(nf);
    at = fft2(b,mf,nf);
    disp(at);
    bt = fft2(a,mf,nf);
  elseif strcmp(pad,'no')
    at = fft2(b);
    bt = fft2(a);
    disp(bt);
  else
    error('Wrong input to XCORRF2');
  end
  
  % % multiply transforms then inverse transform
  % c = ifft2(at.*bt);
  
  % % make real output for real input
  % if ~any(any(imag(a))) && ~any(any(imag(b)))
  %   c = real(c);
  % end
  
  % if strcmp(pad,'yes');    % trim to standard size
  %   c(ma+mb:mf,:) = [];
  %   c(:,na+nb:nf) = [];
  %   elseif strcmp(pad,'no');
  %   c=(c(1:end-1,1:end-1));
  % end
  c = 0;  % Just a filler variable. DELETE
end


% ------ TEST ZONE ------
a = [16 2 3 13 2;
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

[mb, nb] = size(b);
[ma, na] = size(a);
disp("mb: " + mb + newline + "nb: " + nb);
disp("ma: " + ma + newline + "na: " + na);
b = conj(b(mb:-1:1,nb:-1:1));
disp("Reverse Conjugate: ");
disp(b);

mf = 2^nextpow2(ma+mb);
nf = 2^nextpow2(na+nb);
disp("mf: " + mf + newline + "nf: " + nf);

at = fft2(b,mf,nf);
disp("FFT2 of b");
disp(at);

% ------ TEST ZONE ------
