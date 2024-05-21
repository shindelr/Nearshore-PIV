% Called by main, line 36. 
% Leads to firstpass.m (line 20), localfilt.m (line 25),
% naninterp.m (line 26/48), finalpass.m 

function [x,y,u,v,SnR,Pkh]=multipassx(A,B,wins,Dt,overlap,sensit)

A=double(A);
B=double(B);

[sy,sx]=size(A);
iter=size(wins,1);

% Initial passes are for removing large-scale displacements.  Initialize
% displacements (datax,datay) to zero
datax=zeros(floor(sy/(wins(1,1)*(1-overlap))),floor(sx/(wins(1,2)*(1-overlap))));
datay=zeros(floor(sy/(wins(1,1)*(1-overlap))),floor(sx/(wins(1,2)*(1-overlap))));
for i=1:iter-1
  disp(['iter ' num2str(i) ' of ' num2str(iter)])

  % Seems like a good place to use plan_fft in julia, each iteration
  % is a different window size, so do it here. Then again before final
  % pass.

  % PIV
  % func sig: firstpass(A,B,N,overlap,idx,idy)
  [x,y,datax,datay]=firstpass(A,B,wins(i,:),overlap,datax,datay);  %two functions to complete in here.

  % validation.  TODO, simplify these codes!
  % [datax,datay]=globfilt(x,y,datax,datay,3);
  [datax,datay]=localfilt(x,y,datax,datay,sensit,'median',3,[]);
  [datax,datay]=naninterp(datax,datay,'linear',[],x,y);
  datax=floor(datax);
  datay=floor(datay);

  % expand the velocity data to twice the original size
  if(i~=iter-1)
    if wins(i,1)~=wins(i+1,1)
      X=(1:((1-overlap)*2*wins(i+1,1)):sx-2*wins(i+1,1)+1) + wins(i+1,1);
      XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
    else
      XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
      X=XI;
    end
    if wins(i,2)~=wins(i+1,2)
      Y=(1:((1-overlap)*2*wins(i+1,2)):sy-2*wins(i+1,2)+1) + wins(i+1,2);
      YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
    else
      YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
      Y=YI; 
    end
    datax=round(interp2(X,Y',datax,XI,YI'));
    datay=round(interp2(X,Y',datay,XI,YI'));
    [datax,datay]=naninterp(datax,datay,'linear',[],...
                            repmat(XI,size(datax,1),1),...
                            repmat(YI',1,size(datax,2)));  % TODO simplify!
    datax=round(datax);
    datay=round(datay);
  end

end

% Final pass gives displacement to subpixel accuracy
disp('Final iteration')

% Call plan_fft

[x,y,u,v,SnR,Pkh]=finalpass(A,B,wins(end,:),overlap,round(datax),round(datay),Dt);
