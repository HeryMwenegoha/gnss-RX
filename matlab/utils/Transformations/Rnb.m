% Hery A Mwenegoha copyright (c) 2020

function y = Rnb(x)
phi   = x(1);
theta = x(2);
psi   = x(3);
y     = R1phi(phi)*R2theta(theta)*R3psi(psi);
end