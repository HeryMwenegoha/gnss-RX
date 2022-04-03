% Hery A Mwenegoha copyright (c) 2020

function y=Rbw(x)
% Body-frame to wind-frame transformation
% Body-frame: x - along the longitudinal axis
%           : y - along the lateral axis
%           : z - along the normal axis
% wind-frame: x - along the airspeed-vector 
%           : y - laterally, 90-deg orthogonal to the x 
%           : z - completes the right hand orthogonal set
alpha = x(1);
beta  = x(2);
Rbeta   = [ cos(beta),  sin(beta),   0;
           -sin(beta),  cos(beta),   0;
                    0,          0,  1];
RalphaT = [cos(alpha),   0,   sin(alpha);
                    0,   1,            0;
          -sin(alpha),   0,   cos(alpha)];
y = Rbeta * RalphaT;
end