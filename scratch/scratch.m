function [xp,yp,up,vp,SnR,Pkh]=finalpass(A,B,N,ol,idx,idy,Dt)
  if length(N)==1
    M=N;
  else
    M=N(1);
    N=N(2);
  end
  cj=1;
  [sy,sx]=size(A);
  
  % preallocate
  xp=zeros(ceil((size(A,1)-N)/((1-ol)*N))+1, ...
      ceil((size(A,2)-M)/((1-ol)*M))+1);
  yp=xp;
  up=xp;
  vp=xp;
  SnR=xp;
  Pkh=xp;
%   % main loop
  for jj=1:((1-ol)*N):sy-N+1
    ci=1;
    for ii=1:((1-ol)*M):sx-M+1

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

      D2=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));
      E=A(jj:jj+N-1, ii:ii+M-1);
      stad1=std(E(:));
      stad2=std(D2(:));
      if stad1==0
        stad1=1;
      end
      if stad2==0
        stad2=1;
      end
      E=E-mean(E(:));
      F=D2-mean(D2(:));
      %% E(E<0)=0; F(F<0)=0;
  
      % calculate normalized correlation
      R = xcorrf2(E, F) ./ (N * M * stad1 * stad2);
  
      % Find position of maximal value of R
      if all(~isnan(R(:))) && ~all(R(:)==0)  %~isnan(stad1) & ~isnan(stad2)
        if size(R,1)==(N-1)
          [max_y1,max_x1]=find(R==max(R(:)));
        else
          [max_y1,max_x1]=find(R==max(max(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3))));
        end

        if length(max_x1)>1
          max_x1=round(sum(max_x1.^2) ./ sum(max_x1));
          max_y1=round(sum(max_y1.^2)./sum(max_y1));
        end

        if max_x1 == 1 || max_y1 == 1 
          disp([max_x1, max_y1]);
        end

        if max_x1==1
          max_x1=2;
        end
        if max_y1==1
          max_y1=2;
        end
        
      %   % 3-point peak fit using centroid, gaussian (default)
      %   % or parabolic fit
      %   [x0 y0]=intpeak(max_x1,max_y1,R(max_y1,max_x1),...
      %                   R(max_y1,max_x1-1),R(max_y1,max_x1+1),...
      %                   R(max_y1-1,max_x1),R(max_y1+1,max_x1),2,[M,N]);
        
  
      %   % calculate signal to noise ratio
      %   R2=R;
      %   try
      %     %consider changing this from try-catch to a simpler
      %     %distance check. The key here is the distance tot he
      %     %image edge. When peak is close to edge, this NaN
      %     %allocation may fail.
      %     R2(max_y1-3:max_y1+3,max_x1-3:max_x1+3)=NaN;
      %   catch
      %     R2(max_y1-1:max_y1+1,max_x1-1:max_x1+1)=NaN;
      %   end
      %   if size(R,1)==(N-1)
      %     [p2_y2,p2_x2]=find(R2==max(R2(:)));                    
      %   else
      %     [p2_y2,p2_x2]=find(R2==max(max(R2(0.5*N:1.5*N-1,0.5*M:1.5*M-1))));
      %   end
      %   if length(p2_x2)>1
      %     p2_x2=p2_x2(round(length(p2_x2)/2));
      %     p2_y2=p2_y2(round(length(p2_y2)/2));
      %   elseif isempty(p2_x2)
          
      %   end
      %   snr=R(max_y1,max_x1)/R2(p2_y2,p2_x2);
      %   % signal to mean:
      %   %snr=R(max_y1,max_x1)/mean(R(:));
      %   % signal to median:
      %   %snr=R(max_y1,max_x1)/median(median(R(0.5*N+2:1.5*N-3,...
      %   %    0.5*M+2:1.5*M-3)));
  
      %   % store displacements, SnR and Peak Height
      %   up(cj,ci)=(-x0+idx(cj,ci))/Dt;
      %   vp(cj,ci)=(-y0+idy(cj,ci))/Dt;
      %   xp(cj,ci)=(ii+(M/2)-1);
      %   yp(cj,ci)=(jj+(N/2)-1);
      %   SnR(cj,ci)=snr;
      %   Pkh(cj,ci)=R(max_y1,max_x1);
      
      % else
      %   up(cj,ci)=NaN;
      %   vp(cj,ci)=NaN;
      %   SnR(cj,ci)=NaN;
      %   Pkh(cj,ci)=0;
      %   xp(cj,ci)=(ii+(M/2)-1);
      %   yp(cj,ci)=(jj+(N/2)-1);
      end
  
      ci=ci+1;
    end
    cj=cj+1;
  end

%   % now we inline the function XCORRF2 to shave off some time.
%   function c = xcorrf2(a,b,pad)
%   %  c = xcorrf2(a,b)
%   %   Two-dimensional cross-correlation using Fourier transforms.
%   %       XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
%   %       XCORRF2(A) is the autocorrelation function.
%   %       This routine is functionally equivalent to xcorr2 but usually faster.
%   %       See also XCORR2.
  
%   %       Author(s): R. Johnson
%   %       $Revision: 1.0 $  $Date: 1995/11/27 $
  
%   if nargin==2
%       pad='yes';
%   end
  
%   [ma,na] = size(a);
%   %   if nargin == 1
%   %     %       for autocorrelation
%   %     b = a;
%   %   end
%   [mb,nb] = size(b);
%   %       make reverse conjugate of one array
%   b = conj(b(mb:-1:1,nb:-1:1));
%   if strcmp(pad,'yes');
%       %       use power of 2 transform lengths
%       mf = 2^nextpow2(ma+mb);
%       nf = 2^nextpow2(na+nb);
%       at = fft2(b,mf,nf);
%       bt = fft2(a,mf,nf);
%   elseif strcmp(pad,'no');
%       disp("pad no");
%       at = fft2(b);
%       bt = fft2(a);
%   else
%     % disp('Wrong input to XCORRF2'); return
%   end
%   %       multiply transforms then inverse transform
%   c = ifft2(at.*bt);
%   %       make real output for real input
%   if ~any(any(imag(a))) && ~any(any(imag(b)))
%       c = real(c);
%   end
%   if strcmp(pad,'yes');
%       %  trim to standard size
%       c(ma+mb:mf,:) = [];
%       c(:,na+nb:nf) = [];
%   elseif strcmp(pad,'no');
%       c=(c(1:end-1,1:end-1));
      
%       %    c(ma+mb:mf,:) = [];
%       %    c(:,na+nb:nf) = [];
%   end
% end
end
% -----------------------------------------------
% -----------------------------------------------
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
  
      % if cj == 1 & ci == 190
      %   disp(size(idx));
      % end

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
  writematrix(U2, "../tests/mlabOut/first_localfilt/U2.csv");

  teller = 1; 
  [ma, na] = size(U2); 
  histo = zeros(size(nu));

  histostd = zeros(size(nu));
  hista = zeros(size(nu));
  histastd = zeros(size(nu));

  fprintf([' Local ',stat,' filter running: \n'])

  % printed = 0;

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
              % if printed <= 10
              %   disp(["====================="])
              %   disp(["Here on ", ii, jj]);
              %   disp(["usum: ", usum]);
              %   disp(["histostd[j, i]: ", histostd(jj, ii)]);
              %   disp(["====================="])
              %   printed = printed + 1;
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

  %%%%%%%% Locate gridpoints with a higher value than the threshold 

  %%[cy,cx]=find((real(histo)>threshold*real(histostd) | ...
  %%    imag(histo)>threshold*imag(histostd)));

  
  [cy,cx]=find(real(U2) > real(histo) + threshold * real(histostd) |...
      imag(U2) > imag(histo) + threshold * imag(histostd) |...
      real(U2) < real(histo) - threshold * real(histostd) |...
      imag(U2) < imag(histo) - threshold * imag(histostd));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for jj=1:length(cy)
      % uv2(jj)=u(cy(jj),cx(jj)); vv2(jj)=v(cy(jj),cx(jj));
      % xv2(jj)=x(cy(jj),cx(jj)); yv2(jj)=y(cy(jj),cx(jj));
      % Now we asign NotANumber (NaN) to all the points in the matrix that
      % exceeds our threshold.
      nu(cy(jj),cx(jj))=NaN;  nv(cy(jj),cx(jj))=NaN;
  end
  % writematrix(nu, "../tests/mlabOut/mtestNUFILT.csv");
  % writematrix(nv, "../tests/mlabOut/mtestNVFILT.csv");


  rest=length(cy);

  rest2=sum(isnan(u(:)))-sum(prev(:));
  fprintf([num2str(rest),' vectors changed'])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we check for NaN's and interpolate where they exist
  if any(strcmp(varargin,'interp'))
      if any(isnan(u(:)))
          [nu,nv]=naninterp(nu,nv);
      end
  end
  hu=nu(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
  hv=nv(ceil(m/2):end-floor(m/2),ceil(m/2):end-floor(m/2));
  fprintf('.\n')
end
% -----------------------------------------------
% -----------------------------------------------
function [u,v]=naninterp(u,v,varargin)
  usr=1;
  if ~any(strcmp(varargin,'linear')) & ~any(strcmp(varargin,'weighted'))
      met='linear';
  else
      tm=cellfun('isclass',varargin,'char');
      if sum(tm)==3
          disp('Wrong input to naninterp!'); return
      end
      met=varargin(tm(1));
  end
  
  %%%%%%%%%%%%%%%%%%%%%
  if strcmp(met,'weighted')==1 & ~strcmp(met,'linear')==1
      
      if length(varargin)==4
          maske=varargin{2}; 
          if ischar(maske) & ~isempty(maske)
            maske=load(maske); 
            maske=maske.maske; 
          end  
          if isempty(maske)
              [u,v]=naninterp2(u,v); usr=any(isnan(u(:)));
          else
              xx=varargin{3}; yy=varargin{4};
              while usr~=0
                  in2=zeros(size(xx));
                  for i=1:length(maske)
                      in=inpolygon(xx,yy,maske(i).idxw,maske(i).idyw);
                      in2=in2+double(in);                    
                  end
                  in2=logical(in2);
                  % interpolate NaN's using FILLMISS.M
                  u(in2)=0; v(in2)=0;
                  u=fillmiss(u); v=fillmiss(v);
                  usr=any(isnan(u(:)));
                  u(in2)=nan; v(in2)=nan;
              end
              u(in2)=NaN; v(in2)=NaN;
              numm=size(in2,1)*size(in2,2);
          end
      else
          while usr~=0
              numm=sum(isnan(u(:)));
              % interpolate NaN's using FILLMISS.M
              u=fillmiss(u); v=fillmiss(v);
              usr=any(isnan(u(:)));
          end
      end
      fprintf([' -> ',num2str(numm),' NaN''s interpolated\n'])
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(met,'linear')==1 & ~strcmp(met,'weighted')==1
      
    if length(varargin)==4
        maske=varargin{2};
        if ischar(maske) & ~isempty(maske)
          maske=load(maske);
          maske=maske.maske;
        end  

        % Goes here everytime bc x, y, mask are hardcoded: []
        if isempty(maske)
            [u,v]=naninterp2(u,v);
            usr=any(isnan(u(:)));

        else
            xx=varargin{3}; yy=varargin{4};
            maske=rmfield(maske,'msk'); % this is done to avoid the large matrix 
            %being copied into the next function.
            [u,v]=naninterp2(u,v,maske,xx,yy);
        end
    else
      % Do it without masking?
        while usr~=0
            % interpolate NaN's using NANINTERP2.M
            [u,v]=naninterp2(u,v); usr=any(isnan(u(:)));
        end
    end
  else
      disp('Something is VERY wrong with your input'); return  
  end
end
% -----------------------------------------------
% -----------------------------------------------
function [u,v]=naninterp2(u,v,mask,xx,yy)
  % determine Calling m-file:
  [stru,II]=dbstack;
  if length(stru)>2
      test=stru(3).name;
      I2=findstr(test,'multipass');
      if isempty(I2), I2=0; end
  else
      I2=0;
  end 
  if nargin==2
      [py,px]=find(isnan(u)==1);
  else
      py2=[];px2=[]; ipol2=zeros(size(xx));
      for i=1:size(mask,2)
          if I2~=0
              ipol1=inpolygon(xx,yy,mask(i).idx,mask(i).idy);
              ipol2=ipol2+ipol1;
          else
              ipol1=inpolygon(xx,yy,mask(i).idxw,mask(i).idyw);
              ipol2=ipol2+ipol1;
          end 
      end
      [py,px]=find(isnan(u)==1 & ~ipol2 );     
  end
  % counter = 0;
  numm=size(py);
  [dy,dx]=size(u);
  lp=1;
  tel=1;
  % Now sort the NaN's after how many neighbors they have that are
  % physical values. Then we first interpolate those that have 8
  % neighbors, followed by 7, 6, 5, 4, 3, 2 and 1
  % use SORTROWS to sort the numbers
  
  while ~isempty(py)
  
    % check number of neighbors
      for i=1:length(py)
          %correction if vector is on edge of matrix
          corx1=0; corx2=0; cory1=0; cory2=0;
          if py(i)==1
            cory1=1; 
            cory2=0;
          elseif py(i)==dy
            cory1=0;
            cory2=-1;
          end
          if px(i)==1
            corx1=1;
            corx2=0;
          elseif px(i)==dx 
            corx1=0;
            corx2=-1;
          end

          ma = u( py(i)-1+cory1:py(i)+1+cory2, ...
                  px(i)-1+corx1:px(i)+1+corx2 );

          nei(i,1)=sum(~isnan(ma(:)));
          nei(i,2)=px(i);
          nei(i,3)=py(i);
          % writematrix(nei, "../tests/mlabOut/mtest_NEI.csv");
      end
      % writematrix(nei, "../tests/mlabOut/mtest_presortNEI.csv");

      % now sort the rows of NEI to interpolate the vectors with the
      % fewest spurious neighbors.
      nei=flipud(sortrows(nei,1));
      % writematrix(nei, "../tests/mlabOut/mtest_postsortNEI.csv");

      % reconstruct the sorted outlier-vectors.
      % and only interpolate the first 50% of vectors
      ind=find(nei(:,1)>=8);
      while isempty(ind)
          ind=find(nei(:,1) >= 8 - tel);
          tel=tel+1;
      end
      tel=1;
      py=nei(ind,3);
      px=nei(ind,2);

      for j=1:size(py,1)
          corx1=0; corx2=0; cory1=0; cory2=0;
          if py(j)==1
              cory1=1; cory2=0;
          elseif py(j)==dy
              cory1=0; cory2=-1;
          end
          if px(j)==1
              corx1=1; corx2=0;
          elseif px(j)==dx
              corx1=0; corx2=-1;
          end
          tmpu=u(py(j)-1+cory1:py(j)+1+cory2, px(j)-1+corx1:px(j)+1+corx2);
          tmpv=v(py(j)-1+cory1:py(j)+1+cory2, px(j)-1+corx1:px(j)+1+corx2);

          u(py(j),px(j))=mnanmean(tmpu(:));
          v(py(j),px(j))=mnanmean(tmpv(:));

          if lp>numm(1)
             u(py(j),px(j))=0;
             v(py(j),px(j))=0;
          end
      end 

      tt=length(py);
      
      if nargin==2
          [py,px]=find(isnan(u)==1);  
      else
          %in2=zeros(size(xx));
          %for i=1:length(mask)
          %    in=inpolygon(xx,yy,mask(i).idxw,mask(i).idyw);
          %    in2=in2+double(in);
          %end
          [py,px]=find(isnan(u)==1 & ~ipol2 );
      end
      lp=lp+1;

  end
end
% -----------------------------------------------
% -----------------------------------------------


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
Dt = 1; overlap = 0.5; validvec = 3;
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
%     disp(['iter ' num2str(i) ' of ' num2str(iter)]);
    i = 1;
    [x,y,datax,datay] = firstpass(A, B, wins(i, :), overlap, datax, datay);
%     % validation
    [datax,datay]=localfilt(x,y,datax,datay, sensit,'median',3,[]);
    % writematrix(datax, "../tests/mlabOut/multipass_loop/localfilt_datax.csv");
%     [datax,datay]=naninterp(datax,datay,'linear',[],x,y);
    
%     datax=floor(datax);
%     datay=floor(datay);

%     % % expand the velocity data to twice the original size
%     if(i~=iter-1)
%       if wins(i,1)~=wins(i+1,1)
%         X=(1:((1-overlap)*2*wins(i+1,1)):sx-2*wins(i+1,1)+1) + wins(i+1,1);
%         XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
%       else
%         XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
%         X=XI;
%       end
%       if wins(i,2)~=wins(i+1,2)
%         Y=(1:((1-overlap)*2*wins(i+1,2)):sy-2*wins(i+1,2)+1) + wins(i+1,2);
%         YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
%       else
%         YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
%         Y=YI; 
%       end

%       datax=round(interp2(X, Y', datax, XI, YI'));
%       datay=round(interp2(X, Y', datay, XI, YI'));

%       [datax,datay]=naninterp(datax, datay, 'linear', [], ...
%                               repmat(XI, size(datax, 1), 1), ...
%                               repmat(YI', 1, size(datax, 2)) ...
%                               ); 
      
%       datax=round(datax);
%       datay=round(datay);
      

%     end
% end


% % % Final pass gives displacement to subpixel accuracy
% disp('Final iteration')

% % % Call plan_fft

% [x,y,u,v,SnR,Pkh]=finalpass(A,B,wins(end,:),overlap,round(datax),round(datay),Dt);

% ------ TEST ZONE ------
