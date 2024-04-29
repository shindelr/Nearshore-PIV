% Minimal PIV calculation for testing and development.  Uses a sample
% Hopkins 2023 ROXSI image-pair.
%
clear

doprofile=0;  % set to '1' to use the matlab code profiler

% load the two sample image frames into memory
im1=imread('data/im1.jpg');
im2=imread('data/im2.jpg');

% set window sizes.  The algorithm will start by analyzing the frame-pairs
% with an analysis window of size pass_sizes(1), then refines its estimate
% with a smaller window (pass_sizes(2)), and so on.
pivwin=16;  % 16pix is a final window size that worked well for ROXSI 2023
log2pivwin=log2(pivwin);
if(log2pivwin-round(log2pivwin)~=0)
  error('pivwin must be a factor of 2')
end
pass_sizes=2.^[6:-1:log2pivwin]';  % step-down the window sizes in powers of two

% make the analysis window sizes 2d (as required for the matpiv code), and
% repeat the last window iteration (also required)
pass_sizes=[pass_sizes pass_sizes];
pass_sizes=[pass_sizes;pass_sizes(end,:)];

% other input params for matpiv
dt=1;%1/40  % frame time step in seconds.  Set to 1 to get 'pixels per frame' velocity
overlap=0.5;  % fraction of window overlap
validvec=3;  % threshold for vector validation during multipass

% run the PIV analysis
if(doprofile)
  profile on  % optional, initiate the matlab profiler
end
[x,y,u,v,SnR,Pkh]=multipassx(im1,im2,pass_sizes,...
                             dt,overlap,validvec);
if(doprofile)
  profile report  % optional
end

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

% plot the results
clf
pause(.5)
subplot(211)
% pcolor(u),sf      % error: var sf undefined? 
pcolor(u), shading flat
title('u [pixels/frame]')
subplot(212)
% pcolor(v),sf
pcolor(v), shading flat
title('v [pixels/frame]')
for i=1:2
  subplot(2,1,i)
  axis equal tight
  ylim([0 200])
disp('EXITING MAIN')
end
