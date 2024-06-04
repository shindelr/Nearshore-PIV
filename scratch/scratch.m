
function [xx,yy,datax,datay]=firstpass(A,B,N,overlap,idx,idy)

  M=N(1); N=N(2); 
  [sy,sx]=size(A);
  xx=zeros(ceil((size(A,1)-N)/((1-overlap)*N))+1, ...
    ceil((size(A,2)-M)/((1-overlap)*M)) +1);
  yy=xx;
  datax=xx;
  datay=xx; 
  IN=zeros(size(A));

  printed = 0;
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
A=imread('data/im1.jpg');
B=imread('data/im2.jpg');

% JUST SET UP FOR FIRSTPASS/MULTIPASSX
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
  % disp(['iter ' num2str(i) ' of ' num2str(iter)]);
  [x,y,datax,datay] = firstpass(A, B, wins(1, :), overlap, datax, datay);
  writematrix(x, 'mtestX.csv');
  writematrix(y, 'mtestY.csv');
  writematrix(datax, 'mtestDATAX.csv');
  writematrix(datay, 'mtestDATAY.csv');
% end

% ------ TEST ZONE ------
