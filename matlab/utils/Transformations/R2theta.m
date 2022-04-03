% Hery A Mwenegoha copyright (c) 2020

function y = R2theta(theta)
% Rotation around the y-axis: CCW
 y = [ cos(theta),   0,  -sin(theta);
                0,   1,            0;
       sin(theta),   0,   cos(theta)];
end