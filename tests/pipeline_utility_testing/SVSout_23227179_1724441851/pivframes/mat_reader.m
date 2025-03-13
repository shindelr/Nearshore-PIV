
% N = 2
mat355 = load('n2/000355_1724441893930713056.jpg.mat');
mat371 = load('n2/000371_1724441894130751488.jpg.mat');
mat387 = load('n2/000387_1724441894330778968.jpg.mat');

mat355.fn
mat371.fn
mat387.fn

% N = 3
% mat355 = load('n3/000355_1724441893930713056.jpg.mat'); 
% mat372 = load('n3/000372_1724441894143253200.jpg.mat');

% mat355.fn
% mat372.fn


 clf
pause(.5)
subplot(211)
pcolor(mat355.u), shading flat
title('u [pixels/frame]')
subplot(212)
pcolor(mat355.v), shading flat
title('v [pixels/frame]')
for i=1:2
   subplot(2,1,i)
   axis equal tight
   ylim([0 200])
end
