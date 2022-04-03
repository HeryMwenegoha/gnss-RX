% Hery A Mwenegoha copyright (c) 2020

function y = R1phi(phi)
% Rotation around the body x-axis: CW
y  = [1,         0,         0;
      0,  cos(phi),  sin(phi);
      0, -sin(phi),  cos(phi)];
end