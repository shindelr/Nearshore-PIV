max_x1 = 16;
max_y1 = 16;
R = 0.1684079708669066;
Rxm1 =  -8.42042795691657e-17;
% Rxm1 = 0.0000;
Rxp1 = 0.022767921593265242;
Rym1 = 0.013138210493162958;
Ryp1 = -0.020155209279284103;
N = 16.0;
M = 16.0;

% 16.0000   16.0000    0.1684    0.0000    0.0228    0.0131   -0.0202    2.0000   16.0000   16.0000

function [x0,y0]=intpeak(x1,y1,R,Rxm1,Rxp1,Rym1,Ryp1,method,N)
if length(N)==2
    M=N(1); N=N(2);
else
    M=N;
end

if any(find(([R Rxm1 Rxp1 Rym1 Ryp1])==0))
    % to avoid Log of Zero warnings
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


[x_0, y_0] = intpeak(max_x1, max_y1, R, Rxm1, Rxp1, Rym1, Ryp1, 2, [M,N])