
function [xx,yy,datax,datay]=firstpass(A,B,N,overlap,idx,idy)

  M=N(1); N=N(2); 
  [sy,sx]=size(A);
  xx=zeros(ceil((size(A,1)-N)/((1-overlap)*N))+1, ...
      ceil((size(A,2)-M)/((1-overlap)*M)) +1);
  yy=xx;
  datax=xx;
  datay=xx; 
  IN=zeros(size(A));

  
  cj=1;
  for jj=1:((1-overlap)*N):sy-N+1
    ci=1;
    for ii=1:((1-overlap)*M):sx-M+1 
  
      if IN(jj+N/2,ii+M/2)~=1 
        
        if isnan(idx(cj,ci))
          idx(cj,ci)=0;
        end
        if isnan(idy(cj,ci))
          idy(cj,ci)=0;
        end
        if jj+idy(cj,ci)<1
          idy(cj,ci)=1-jj;
        elseif jj+idy(cj,ci)>sy-N+1
          idy(cj,ci)=sy-N+1-jj;
        end       
        if ii+idx(cj,ci)<1
          idx(cj,ci)=1-ii;    
        elseif ii+idx(cj,ci)>sx-M+1
          idx(cj,ci)=sx-M+1-ii;
        end
  
        C=A(jj:jj+N-1,ii:ii+M-1);   
        D=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));
        C=C-mean(C(:)); D=D-mean(D(:)); %C(C<0)=0; D(D<0)=0;
        stad1=std(C(:)); stad2=std(D(:)); 
        %C=C.*W; %D=D.*W;  % Apply weight function by uncommenting
  
        if stad1==0, stad1=nan; end
        if stad2==0, stad2=nan; end
        
        % normalized correlation
        R=xcorrf2(C,D)/(N*M*stad1*stad2);
  
        % find position of maximal value of R
        if size(R,1)==(N-1)
          [max_y1,max_x1]=find(R==max(R(:)));
        else
          [max_y1,max_x1]=find(R==max(max(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3))));
        end
        disp(max_y1);
        disp(max_x1);
        
        if length(max_x1)>1
          max_x1=round(sum(max_x1.*(1:length(max_x1))')./sum(max_x1));
          max_y1=round(sum(max_y1.*(1:length(max_y1))')./sum(max_y1));
        elseif isempty(max_x1)
          idx(cj,ci)=nan;
          idy(cj,ci)=nan;
          max_x1=nan;
          max_y1=nan;
        end

        % store the displacements in variable datax/datay
        datax(cj,ci)=-(max_x1-(M))+idx(cj,ci);
        datay(cj,ci)=-(max_y1-(N))+idy(cj,ci);
        xx(cj,ci)=ii+M/2; yy(cj,ci)=jj+N/2;
        ci=ci+1;
      else
        xx(cj,ci)=ii+M/2; yy(cj,ci)=jj+N/2;
        datax(cj,ci)=NaN; datay(cj,ci)=NaN; ci=ci+1;
      end
    end
  
    cj=cj+1;
  end
end

% -----------------------------------------------
% -----------------------------------------------
function c = xcorrf2(a,b,pad)
  if nargin==2
    pad='yes';
  end

  [ma,na] = size(a);
  [mb,nb] = size(b);

  % make reverse conjugate of one array
  b = conj(b(mb:-1:1,nb:-1:1));
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

  % multiply transforms then inverse transform
  c = ifft2(at.*bt);

  % make real output for real input
  if ~any(any(imag(a))) && ~any(any(imag(b)))
    c = real(c);
  end
  if strcmp(pad,'yes');    % trim to standard size
    c(ma+mb:mf,:) = [];
    c(:,na+nb:nf) = [];
  elseif strcmp(pad,'no');
    c=(c(1:end-1,1:end-1));
  end
end


% ------ TEST ZONE ------
% A = [
%     38, 28, 14, 42, 7, 20, 38, 18, 22, 10;
%     10, 23, 35, 39, 23, 2, 21, 1, 23, 43;
%     29, 37, 1, 20, 32, 11, 21, 43, 24, 48;
%     26, 41, 27, 15, 14, 46, 50, 43, 2, 36;
%     50, 6, 20, 8, 38, 17, 3, 24, 13, 49;
%     8, 25, 1, 19, 27, 46, 6, 43, 7, 46;
%     34, 13, 16, 35, 49, 39, 3, 1, 5, 41;
%     3, 28, 17, 25, 43, 33, 9, 35, 13, 30;
%     47, 14, 7, 13, 22, 39, 20, 15, 44, 17;
%     46, 23, 25, 24, 44, 40, 28, 14, 44, 0
% ];

% B = [
%     46, 34, 22, 49, 7, 27, 45, 28, 24, 10;
%     17, 25, 37, 39, 33, 6, 30, 7, 32, 51;
%     35, 45, 8, 21, 32, 17, 27, 50, 28, 50;
%     33, 46, 37, 17, 14, 48, 54, 45, 2, 40;
%     59, 12, 26, 18, 46, 26, 12, 26, 19, 49;
%     11, 28, 5, 25, 33, 56, 9, 49, 17, 48;
%     39, 14, 25, 43, 53, 44, 6, 11, 14, 47;
%     11, 34, 17, 25, 51, 43, 17, 38, 21, 32;
%     53, 19, 14, 23, 30, 43, 20, 17, 53, 24;
%     56, 28, 32, 32, 47, 40, 28, 23, 47, 6
% ];

% A = rand(2048, 3072);
% B = A + 2;
A=imread('data/im1.jpg');
B=imread('data/im2.jpg');

% JUST SET UP FOR FIRSTPASS
% ------------------------------------------------------------------------------
pivwin=16;
log2pivwin=log2(pivwin);
pass_sizes=2.^[6:-1:log2pivwin]';
pass_sizes=[pass_sizes pass_sizes];
pass_sizes=[pass_sizes;pass_sizes(end,:)];
wins = pass_sizes; 
dt = 1; overlap = 0.5; validvec = 3;
A=double(A);
B=double(B);
[sy,sx]=size(A);
iter=size(wins,1);
datax=zeros(floor(sy/(wins(1,1)*(1-overlap))),floor(sx/(wins(1,2)*(1-overlap))));
datay=zeros(floor(sy/(wins(1,1)*(1-overlap))),floor(sx/(wins(1,2)*(1-overlap))));
% ------------------------------------------------------------------------------

% for i=1:iter-1
%   disp(['iter ' num2str(i) ' of ' num2str(iter)]);
x = firstpass(A, B, wins(1, :), overlap, datax, datay);
%   disp('Done');
% end

% ------ TEST ZONE ------
