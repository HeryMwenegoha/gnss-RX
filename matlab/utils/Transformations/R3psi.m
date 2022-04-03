% Hery A Mwenegoha copyright (c) 2020

function y = R3psi(psi)
% Rotation around the z-axis: CW
y = [ cos(psi),    sin(psi),    0;
     -sin(psi),    cos(psi),    0;
             0,           0,    1];
end