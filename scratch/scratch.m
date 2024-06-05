
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
% -----------------------------------------------
% -----------------------------------------------
function [hu,hv]=localfilt(x,y,u,v,threshold,varargin)
  IN=zeros(size(u));

  if nargin < 5
      disp(' Not enough input arguments!'); return
  end
  if ischar(threshold)
      disp(' Please give threshold as numeric input'); return
  end

  if nargin > 5
      tm=cellfun('isclass',varargin,'double'); 
      pa=find(tm==1); 
      if ~isempty(pa)
        m=cat(1,varargin{pa}); 
      else
        m=3; 
      end 
      if any(strcmp(varargin,'median'))
          method='mnanmedian'; stat='median'; ff=1;

      elseif any(strcmp(varargin,'mean'))
          method='mnanmean'; stat='mean'; ff=2;
      end  
      if ~any(strcmp(varargin,'mean')) & ~any(strcmp(varargin,'median'))
          method='mnanmedian'; stat='median';
      end
      if nargin==8
          maske=varargin{end}; 
          if ischar(maske) & ~isempty(maske)
            maske=load(maske); 
            maske=maske.maske; 
          end

          if ~isempty(maske)
            for ii=1:length(maske) 
              IN2=inpolygon(x,y,maske(ii).idxw,maske(ii).idyw);
              IN=[IN+IN2];
            end
            
          else 
            IN=zeros(size(u)); 
          end
      end
  end

  if nargin==5
      m=3; method='mnanmedian'; stat='median';
  end
  nu = zeros(size(u) + 2 * floor(m/2)) * nan;
  nv = zeros(size(u) + 2 * floor(m/2)) * nan;
  nu(floor(m/2) + 1:end - floor(m/2), floor(m/2) + 1:end - floor(m/2)) = u;
  nv(floor(m/2) + 1:end - floor(m/2), floor(m/2) + 1:end - floor(m/2)) = v;
  
  INx=zeros(size(nu));
  INx(floor(m/2) + 1:end - floor(m/2), floor(m/2) + 1:end - floor(m/2)) = IN;
  
  prev = isnan(nu); 
  previndx = find(prev==1); 
  U2 = nu + i * nv;
  % writematrix(U2, "../tests/mlabOut/mtestU2.csv");  

  teller = 1; 
  [ma, na] = size(U2); 
  histo = zeros(size(nu));

  histostd = zeros(size(nu));
  hista = zeros(size(nu));
  histastd = zeros(size(nu));

  fprintf([' Local ',stat,' filter running: \n'])

  printed = 0;

  for ii=m-1:1:na-m+2  
      for jj=m-1:1:ma-m+2
          if INx(jj,ii)~=1
              
              tmp=U2(jj-floor(m/2):jj+floor(m/2), ii-floor(m/2):ii+floor(m/2)); 
              tmp(ceil(m/2),ceil(m/2))=NaN;

              
              if ff==1
                usum=mnanmedian(tmp(:));
              elseif ff==2
                usum=mnanmean(tmp(:));
              end
              histostd(jj,ii)=mnanstd(tmp(:));
              % if printed == 0
              %   disp(tmp(:));
              %   disp(histostd(jj, ii));
              %   printed = printed + 1;
              %   % disp(size(tmp));
              % end
          else
              usum=nan; tmp=NaN; histostd(jj,ii)=nan;
          end
          % u1=real(usum).^2 - real(U2(jj,ii)).^2;
          % v1=imag(usum).^2 - imag(U2(jj,ii)).^2;
          
          % histo(jj,ii)=u1+i*v1;
          histo(jj,ii)=usum;
          % histostd(jj,ii)=mnanstd(real(tmp(:))) + i*mnanstd(imag(tmp(:)));
          
          % th1=angle(usum); th2=angle(U2(jj,ii));
          % if th1<0, th1=2*pi+th1; end
          % if th2<0, th2=2*pi+th2; end
          % hista(jj,ii)=(th1-th2);
          % if hista(jj,ii)<0, hista(jj,ii)=2*pi+hista(jj,ii); end 
          % histastd(jj,ii)=mnanstd(abs(angle(tmp(:))));
      end
      fprintf('.')
      
  end

  % %%%%%%%% Locate gridpoints with a higher value than the threshold 

  % %[cy,cx]=find((real(histo)>threshold*real(histostd) | ...
  % %    imag(histo)>threshold*imag(histostd)));
  % [cy,cx]=find( ( real(U2)>real(histo)+threshold*real(histostd) |...
  %     imag(U2)>imag(histo)+threshold*imag(histostd) |...
  %     real(U2)<real(histo)-threshold*real(histostd) |...
  %     imag(U2)<imag(histo)-threshold*imag(histostd) ) );

  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % for jj=1:length(cy)
  %     %uv2(jj)=u(cy(jj),cx(jj)); vv2(jj)=v(cy(jj),cx(jj));
  %     %xv2(jj)=x(cy(jj),cx(jj)); yv2(jj)=y(cy(jj),cx(jj));
  %     % Now we asign NotANumber (NaN) to all the points in the matrix that
  %     % exceeds our threshold.
  %     nu(cy(jj),cx(jj))=NaN;  nv(cy(jj),cx(jj))=NaN;
  % end

  % rest=length(cy);

  % rest2=sum(isnan(u(:)))-sum(prev(:));
  % fprintf([num2str(rest),' vectors changed'])
  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % % Now we check for NaN's and interpolate where they exist
  % if any(strcmp(varargin,'interp'))
  %     if any(isnan(u(:)))
  %         [nu,nv]=naninterp(nu,nv);
  %     end
  % end
  % hu=nu(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
  % hv=nv(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
  % fprintf('.\n')

  % Dummies
  hu = 0;
  hv = 0;
end

% ------ TEST ZONE ------
tic
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
sensit = 3;
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% function Multipassx technically
% Loop disabled for testing
% for i=1:iter-1
    % disp(['iter ' num2str(i) ' of ' num2str(iter)]);
    [x,y,datax,datay] = firstpass(A, B, wins(1, :), overlap, datax, datay);

  % validation
  [datax,datay]=localfilt(x,y,datax,datay, sensit,'median',3,[]);
  toc
  % [datax,datay]=naninterp(datax,datay,'linear',[],x,y);
  % datax=floor(datax);
  % datay=floor(datay);

  % % expand the velocity data to twice the original size
  % if(i~=iter-1)
  %   if wins(i,1)~=wins(i+1,1)
  %     X=(1:((1-overlap)*2*wins(i+1,1)):sx-2*wins(i+1,1)+1) + wins(i+1,1);
  %     XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
  %   else
  %     XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
  %     X=XI;
  %   end
  %   if wins(i,2)~=wins(i+1,2)
  %     Y=(1:((1-overlap)*2*wins(i+1,2)):sy-2*wins(i+1,2)+1) + wins(i+1,2);
  %     YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
  %   else
  %     YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
  %     Y=YI; 
  %   end
  %   datax=round(interp2(X,Y',datax,XI,YI'));
  %   datay=round(interp2(X,Y',datay,XI,YI'));
  %   [datax,datay]=naninterp(datax,datay,'linear',[],...
  %                           repmat(XI,size(datax,1),1),...
  %                           repmat(YI',1,size(datax,2)));  % TODO simplify!
  %   datax=round(datax);
  %   datay=round(datay);
  % end
  % end
  % writematrix(x, 'mtestX.csv');
  % writematrix(y, 'mtestY.csv');
  % writematrix(datax, "../tests/mlabOut/mtestDATAX.csv");
  % writematrix(datay, "../tests/mlabOut/mtestDATAY.csv");

% ------ TEST ZONE ------
