function c = xcorrf2(a,b,pad)
if nargin==2
    pad='yes';
c = 4;

%   [ma,na] = size(a);
%   [mb,nb] = size(b);

%   % make reverse conjugate of one array
%   b = conj(b(mb:-1:1,nb:-1:1));
%   if strcmp(pad,'yes');   % use power of 2 transform lengths
%     mf = 2^nextpow2(ma+mb);
%     nf = 2^nextpow2(na+nb);
%     at = fft2(b,mf,nf);
%     bt = fft2(a,mf,nf);
%   elseif strcmp(pad,'no');
%     at = fft2(b);
%     bt = fft2(a);
%   else
%     error('Wrong input to XCORRF2');
%   end

%   % multiply transforms then inverse transform
%   c = ifft2(at.*bt);

%   % make real output for real input
%   if ~any(any(imag(a))) && ~any(any(imag(b)))
%     c = real(c);
%   end

%   if strcmp(pad,'yes');    % trim to standard size
%     c(ma+mb:mf,:) = [];
%     c(:,na+nb:nf) = [];
%     elseif strcmp(pad,'no');
%     c=(c(1:end-1,1:end-1));
%   end
end


% ------ TEST ZONE ------
c = xcorrf2(1);

% ------ TEST ZONE ------
