% Hery A Mwenegoha copyright (c) 2020

function y = euler2quat(x)
% Euler to Quaternion
phi   = x(1);
theta = x(2);
psi   = x(3);
q0    = sin(phi/2).*sin(psi/2).*sin(theta/2) +...
    cos(phi/2).*cos(psi/2).*cos(theta/2);
q1    = cos(psi/2).*cos(theta/2).*sin(phi/2) -...
    cos(phi/2).*sin(psi/2).*sin(theta/2);
q2    = cos(phi/2).*cos(psi/2).*sin(theta/2) +...
    cos(theta/2).*sin(phi/2).*sin(psi/2);
q3    = cos(phi/2).*cos(theta/2).*sin(psi/2) -...
    cos(psi/2).*sin(phi/2).*sin(theta/2);
y     = [q0;q1;q2;q3];
end