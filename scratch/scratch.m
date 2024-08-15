function [hu,hv]=globfilt(x,y,u,v,varargin)
if nargin < 5
  disp('Not enough input arguments!'); return
end
tm=cellfun('isclass',varargin,'double');
pa=find(tm==1);
if length(pa)>1
  disp('Only one numeric input allowed!'); return
end
fprintf(' Global filter running - ')
if max(sqrt(u(:).^2+v(:).^2))~=0
  scale=2/max(sqrt(u(:).^2+v(:).^2));
else
  scale=0.1;
end
if any(strcmp(varargin,'manual')) & ~any(strcmp(varargin,'loop'))
  figure, subplot(211),vekplot2(x,y,u,v,scale);
  subplot(212),plot(u,v,'.'), title('scatter plot of velocities')
  xlabel('worldcoordinates per second')
  ylabel('worldcoordinates per second')
  disp('Use left button to mark the 4 corners around your region...')
  hold on
  for i=1:4
    [ii(i),jj(i)]=ginput(1);
    if i>1
      h1=plot([ii(i-1) ii(i)],[jj(i-1) jj(i)],'k-');
      set(h1,'LineWidth',[2]);
    end
  end
  h2=plot([ii(4) ii(1)],[jj(4) jj(1)],'k-');
  set(h2,'LineWidth',[2]);
  clear in; in=inpolygon(u,v,ii,jj);
  subplot(211), hold on
  vekplot2(x,y,u,v,scale,'b');
  vekplot2(x(~in),y(~in),u(~in),v(~in),scale,'r');
  drawnow
elseif any(strcmp(varargin,'loop')) & ~any(strcmp(varargin,'manual'))
  usr=1;
  if ~isempty(pa)
    param=cat(1,varargin{pa});
    %param=cell2mat(varargin(pa));
  else
    param=3;
    disp('Warning! no threshold specified. Using standard setting.')
  end
  xo=mnanmean(u(:)); yo=mnanmean(v(:));
  while param~=0,
    sx=param*mnanstd(u(:)); sy=param*mnanstd(v(:));
    if ~any(strcmp(varargin,'circle'))
      ii=[xo+sx; xo+sx; xo-sx; xo-sx];
      jj=[yo+sy; yo-sy; yo-sy; yo+sy];
    else
      ttt=0:0.1:2*pi;
      ii=xo+sx*sin(ttt); jj=yo+sy*cos(ttt); ii=ii(:); jj=jj(:);
    end
    figure(gcf), subplot(211),vekplot2(x,y,u,v,scale);
    subplot(212), plot(u,v,'.'), hold on, plot(ii,jj,'-')
    plot([ii(1) ii(4)],[jj(4) jj(4)],'-'),  hold off
    clear in; in=inpolygon(u,v,ii,jj);
    subplot(211), hold on
    vekplot2(x,y,u,v,scale,'b');
    vekplot2(x(~in),y(~in),u(~in),v(~in),scale,'r');
    fprintf(['with limit: ',num2str(param),...
      ' *std [U V]'])
    param=input(['To change THRESHOLD type new value, \n type 0 to use current value >> ']);
    
  end
  close
elseif any(cellfun('isclass',varargin,'double')==1)& ~any(strcmp(varargin,'loop'))
  param=cat(1,varargin{pa});
  %param=cell2mat(varargin(pa));
  if length(param)==1
    sx=param*mnanstd(u(:));
    sy=param*mnanstd(v(:));
    
    xo=mnanmean(u(:));
    yo=mnanmean(v(:));
    
    if ~any(strcmp(varargin,'circle'))
      ii=[xo+sx; xo+sx; xo-sx; xo-sx];
      jj=[yo+sy; yo-sy; yo-sy; yo+sy];
      % else
      %   ttt=0:0.1:2*pi;
      %   ii=xo+sx*sin(ttt); jj=yo+sy*cos(ttt); ii=ii(:); jj=jj(:);
    end
    fprintf(['with limit: ',num2str(param),...
      ' *std [U V]'])
    %close
    % elseif length(param)==4
    %   sx(1)=param(1); sy(1)=param(3);
    %   sx(2)=param(2); sy(2)=param(4);
    %   ii=[sx(1); sx(2); sx(2); sx(1)];
    %   jj=[sy(1); sy(1); sy(2); sy(2)];
    %   fprintf(['Current limit: [',num2str(param),'] = (umin,umax,vmin,vmax)'])
    % else
    %   fprintf('Something wrong with your numerical input')
  end
  % else
  %   disp('Error! Check your input to GLOBFILT')
  %   %close
  %   return
end
prev=isnan(u); previndx=find(prev==1);
%Locate points inside chosen area
in=inpolygon(u,v,ii,jj);
disp(in)
nx=x(~in); ny=y(~in);
nu=u(~in); nv=v(~in);

% %scale=3/max(sqrt(u(:).^2+v(:).^2));
% if any(strcmp(varargin,'manual')==1) | any(strcmp(varargin,'loop')==1)
%   figure, vekplot2(x(:).',y(:).',u(:).',v(:).',scale,'b');
%   hold on, grid on
%   vekplot2(nx(:).',ny(:).',nu(:).',nv(:).',scale,'r');
%   xlabel([num2str(length(nx(:))-length(previndx(:))),...
%     ' outliers identified by this filter, from totally ',...
%     num2str(length(u(:))),' vectors'])
% end
% %Exclude points outside area
% u(~in)=NaN; v(~in)=NaN;
% %interpolate
% if any(strcmp(varargin,'interp')==1)
%   if any(isnan(u(:)))
%     [u,v]=naninterp2(u,v);
%     vekplot2(x,y,u,v,scale,'g');
%     title('Green arrows are validated and interpolated vector field')
%   end
% end
% fprintf([' ..... ',num2str(length(nx(:))-length(previndx(:))),...
%       ' vectors changed\n'])
% hu=u; hv=v;
end
% -----------------------------------------------
% -----------------------------------------------
function [x0,y0]=intpeak(x1,y1,R,Rxm1,Rxp1,Rym1,Ryp1,method,N)
if length(N)==2
  M=N(1); N=N(2);
else
  M=N;
end

if any(find(([R Rxm1 Rxp1 Rym1 Ryp1])==0))
  % to avoid Log of Zero warnings
  disp("Here")
  method=1;
end

if method==1
  x01=(((x1-1)*Rxm1)+(x1*R)+((x1+1)*Rxp1)) / (Rxm1+ R+Rxp1);
  y01=(((y1-1)*Rym1)+(y1*R)+((y1+1)*Ryp1)) / (Rym1+ R+Ryp1);
  x0=x01-(M);
  y0=y01-(N);
elseif method==2
  x01=x1 + ( (log(Rxm1)-log(Rxp1))/( (2*log(Rxm1))-(4*log(R))+(2*log(Rxp1))) );
  y01=y1 + ( (log(Rym1)-log(Ryp1))/( (2*log(Rym1))-(4*log(R))+(2*log(Ryp1))) );
  x0=x01-(M);
  y0=y01-(N);
elseif method==3
  x01=x1 + ( (Rxm1-Rxp1)/( (2*Rxm1)-(4*R)+(2*Rxp1)) );
  y01=y1 + ( (Rym1-Ryp1)/( (2*Rym1)-(4*R)+(2*Ryp1)) );
  x0=x01-(M);
  y0=y01-(N);
  
  
else
  
  disp(['Please include your desired peakfitting function; 1 for',...
    ' 3-point fit, 2 for gaussian fit, 3 for parabolic fit'])
  
end


x0=real(x0);
y0=real(y0);
end
% -----------------------------------------------
% -----------------------------------------------
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
% main loop
% fid = fopen("../tests/mlabOut/finalpass/intpeak.csv", "a");
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
    
    % if cj == 196 && ci == 46
    %   disp([max_y1, max_x1])
    %   writematrix(R, "../tests/mlabOut/finalpass/R_cj196_ci46.csv");
    % end
    
    % Find position of maximal value of R
    if all(~isnan(R(:))) && ~all(R(:)==0)  %~isnan(stad1) & ~isnan(stad2)
      if size(R,1)==(N-1)
        [max_y1,max_x1]=find(R==max(R(:)));
      else
        [max_y1,max_x1]=find(R==max(max(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3))));
        
        % if cj == 196 && ci == 46
        %   disp([max_y1, max_x1])
        %   writematrix(R(0.5*N+2:1.5*N-3,0.5*M+2:1.5*M-3), "../tests/mlabOut/finalpass/subset.csv");
        % end
        
      end
      
      if length(max_x1)>1
        max_x1=round(sum(max_x1.^2) ./ sum(max_x1));
        max_y1=round(sum(max_y1.^2)./sum(max_y1));
      end
      
      % if max_x1 == 1 || max_y1 == 1
      %   disp([max_x1, max_y1]);
      % end
      
      if max_x1==1
        max_x1=2;
      end
      if max_y1==1
        max_y1=2;
      end
      
      
      % 3-point peak fit using centroid, gaussian (default)
      % or parabolic fit
      [x0 y0]=intpeak(max_x1,max_y1,R(max_y1,max_x1),...
        R(max_y1,max_x1-1),R(max_y1,max_x1+1),...
        R(max_y1-1,max_x1),R(max_y1+1,max_x1),2,[M,N]);
      
      if cj == 196 && ci == 46
        disp([max_y1, max_x1])
        % disp([max_x1,max_y1,R(max_y1,max_x1),...
        %   R(max_y1,max_x1-1),R(max_y1,max_x1+1),...
        %   R(max_y1-1,max_x1),R(max_y1+1,max_x1),2,[M,N]])
      end
      
      % fprintf(fid, "%d, %d, %f, %f\n", cj, ci, x0, y0);
      
      
      % calculate signal to noise ratio
      R2=R;
      
      try
        %consider changing this from try-catch to a simpler
        %distance check. The key here is the distance tot he
        %image edge. When peak is close to edge, this NaN
        %allocation may fail.
        R2(max_y1-3:max_y1+3,max_x1-3:max_x1+3)=NaN;
      catch
        R2(max_y1-1:max_y1+1,max_x1-1:max_x1+1)=NaN;
      end
      
      
      if size(R,1)==(N-1)
        [p2_y2,p2_x2]=find(R2==max(R2(:)));
      else
        [p2_y2,p2_x2]=find(R2==max(max(R2(0.5*N:1.5*N-1,0.5*M:1.5*M-1))));
        % if cj == 1 && ci == 1
        %   disp([p2_y2, p2_x2])
        %   % writematrix(R2(0.5*N:1.5*N-1,0.5*M:1.5*M-1), "../tests/mlabOut/finalpass/subset.csv")
        % end
      end
      
      if length(p2_x2)>1
        p2_x2=p2_x2(round(length(p2_x2)/2));
        p2_y2=p2_y2(round(length(p2_y2)/2));
      elseif isempty(p2_x2)
        disp("Empty set found")
      end
      
      snr=R(max_y1,max_x1)/R2(p2_y2,p2_x2);
      % signal to mean:
      % snr=R(max_y1,max_x1)/mean(R(:));
      % signal to median:
      % snr=R(max_y1,max_x1)/median(median(R(0.5*N+2:1.5*N-3,...
      %    0.5*M+2:1.5*M-3)));
      
      % store displacements, SnR and Peak Height
      up(cj,ci)=(-x0+idx(cj,ci))/Dt;
      % if cj == 1 && ci == 1
      %   % disp([x0, idx(cj,ci), Dt])
      %   disp((x0 + idx(cj,ci)) / Dt)
      % end
      
      vp(cj,ci)=(-y0+idy(cj,ci))/Dt;
      xp(cj,ci)=(ii+(M/2)-1);
      yp(cj,ci)=(jj+(N/2)-1);
      SnR(cj,ci)=snr;
      Pkh(cj,ci)=R(max_y1,max_x1);
      
    else
      up(cj,ci)=NaN;
      vp(cj,ci)=NaN;
      SnR(cj,ci)=NaN;
      Pkh(cj,ci)=0;
      xp(cj,ci)=(ii+(M/2)-1);
      yp(cj,ci)=(jj+(N/2)-1);
    end
    
    ci=ci+1;
  end
  cj=cj+1;
end
% fclose(fid);

% now we inline the function XCORRF2 to shave off some time.
  function c = xcorrf2(a,b,pad)
    %  c = xcorrf2(a,b)
    %   Two-dimensional cross-correlation using Fourier transforms.
    %       XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
    %       XCORRF2(A) is the autocorrelation function.
    %       This routine is functionally equivalent to xcorr2 but usually faster.
    %       See also XCORR2.
    
    %       Author(s): R. Johnson
    %       $Revision: 1.0 $  $Date: 1995/11/27 $
    
    if nargin==2
      pad='yes';
    end
    
    [ma,na] = size(a);
    %   if nargin == 1
    %     %       for autocorrelation
    %     b = a;
    %   end
    [mb,nb] = size(b);
    %       make reverse conjugate of one array
    b = conj(b(mb:-1:1,nb:-1:1));
    if strcmp(pad,'yes');
      %       use power of 2 transform lengths
      mf = 2^nextpow2(ma+mb);
      nf = 2^nextpow2(na+nb);
      at = fft2(b,mf,nf);
      bt = fft2(a,mf,nf);
    elseif strcmp(pad,'no');
      disp("pad no");
      at = fft2(b);
      bt = fft2(a);
    else
      % disp('Wrong input to XCORRF2'); return
    end
    %       multiply transforms then inverse transform
    c = ifft2(at.*bt);
    %       make real output for real input
    if ~any(any(imag(a))) && ~any(any(imag(b)))
      c = real(c);
    end
    if strcmp(pad,'yes');
      %  trim to standard size
      c(ma+mb:mf,:) = [];
      c(:,na+nb:nf) = [];
    elseif strcmp(pad,'no');
      c=(c(1:end-1,1:end-1));
      
      %    c(ma+mb:mf,:) = [];
      %    c(:,na+nb:nf) = [];
    end
  end
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

teller = 1;
[ma, na] = size(U2);
histo = zeros(size(nu));

histostd = zeros(size(nu));
hista = zeros(size(nu));
histastd = zeros(size(nu));

fprintf([' Local ',stat,' filter running: \n'])

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
    else
      usum=nan; tmp=NaN; histostd(jj,ii)=nan;
    end
    
    % u1=real(usum).^2 - real(U2(jj,ii)).^2;
    % v1=imag(usum).^2 - imag(U2(jj,ii)).^2;
    
    % histo(jj,ii)=u1+i*v1;
    
    
    histo(jj,ii)=usum;
    
    % if jj == 5 && ii == 10
    %   disp(["Iter: ", jj, ii])
    %   histo(jj,ii)=usum;
    %   disp(["usum = ", usum])
    % end
    
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
% writematrix(histo, "../tests/mlabOut/first_localfilt/histo.csv");


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
% writematrix(nu, "../tests/mlabOut/first_localfilt/nu.csv");
% writematrix(nv, "../tests/mlabOut/first_localfilt/nv.csv");

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

% writematrix(hu, "../tests/mlabOut/first_localfilt/hu.csv");
% writematrix(hv, "../tests/mlabOut/first_localfilt/hv.csv");

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
% tic;
% for i=1:iter-1
%   % i = 1;
%   disp(['iter ' num2str(i) ' of ' num2str(iter)]);
%   [x,y,datax,datay] = firstpass(A, B, wins(i, :), overlap, datax, datay);
%   % filename = sprintf('firstpass_datax%d.csv', i);
%   % path = sprintf("../tests/mlabOut/multipass_loop/%s", filename);
%   % writematrix(datax, path);

%   % validation
%   [datax,datay]=localfilt(x,y,datax,datay, sensit,'median',3,[]);
%   % writematrix(datax, "../tests/mlabOut/multipass_loop/localfilt_datax.csv");
%   % writematrix(datay, "../tests/mlabOut/multipass_loop/localfilt_datay.csv");

%   [datax,datay]=naninterp(datax,datay,'linear',[],x,y);

%   datax=floor(datax);
%   datay=floor(datay);
%   % writematrix(datax, "../tests/mlabOut/multipass_loop/1stpass_linnaninterp_datax.csv");
%   % writematrix(datay, "../tests/mlabOut/multipass_loop/1stpass_linnaninterp_datay.csv");

%   % expand the velocity data to twice the original size
%   if(i~=iter-1)
%     if wins(i,1)~=wins(i+1,1)
%       X=(1:((1-overlap)*2*wins(i+1,1)):sx-2*wins(i+1,1)+1) + wins(i+1,1);
%       XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
%     else
%       XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
%       X=XI;
%     end
%     if wins(i,2)~=wins(i+1,2)
%       Y=(1:((1-overlap)*2*wins(i+1,2)):sy-2*wins(i+1,2)+1) + wins(i+1,2);
%       YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
%     else
%       YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
%       Y=YI;
%     end

%     datax=round(interp2(X, Y', datax, XI, YI'));
%     datay=round(interp2(X, Y', datay, XI, YI'));

%     [datax,datay]=naninterp(datax, datay, 'linear', [], ...
%       repmat(XI, size(datax, 1), 1), ...
%       repmat(YI', 1, size(datax, 2)) ...
%       );

%     datax=round(datax);
%     datay=round(datay);

%     % writematrix(datax, "../tests/mlabOut/multipass_loop/reginterp_datax.csv");
%     % writematrix(datay, "../tests/mlabOut/multipass_loop/reginterp_datay.csv");
%   end

% end
% writematrix(datax, "../tests/mlabOut/penultimate_datax.csv");
% writematrix(datay, "../tests/mlabOut/penultimate_datay.csv");
%
% % % Final pass gives displacement to subpixel accuracy
disp('Final iteration')
datax = readmatrix("/tests/mlabOut/penultimate_datax.csv");
datay = readmatrix("/tests/mlabOut/penultimate_datay.csv");

[x,y,u,v,SnR,Pkh]=finalpass(A,B,wins(end,:),overlap,round(datax),round(datay),Dt);
% writematrix(x, "../tests/mlabOut/finalpass/x.csv");
% writematrix(y, "../tests/mlabOut/finalpass/y.csv");
% writematrix(u, "../tests/mlabOut/finalpass/u.csv");
% writematrix(v, "../tests/mlabOut/finalpass/v.csv");
% writematrix(SnR, "../tests/mlabOut/finalpass/SnR.csv");
% writematrix(Pkh, "../tests/mlabOut/finalpass/Pkh.csv");
% toc

% save original results prior to quality control...
uraw=u;
vraw=v;

% snrfilt: reject data with too-low signal to noise level
snrthresh=1.3;
ibad=find(SnR<snrthresh);
u(ibad)=nan;
v(ibad)=nan;

% peakfilt: reject data for which the correlation peak was too low
pkhthresh=0.3;
ibad=find(Pkh<pkhthresh);
u(ibad)=nan;
v(ibad)=nan;

% globfilt: reject data that disagree strongly with its neighbors in a local
% window
[u,v]=globfilt(x,y,u,v,3);


% ------ TEST ZONE ------
