% Hery A Mwenegoha copyright (c) 2020

function y = Rnb(varargin)
if nargin == 1 || nargin == 3
    x = [varargin{:}];
else
    error('input is either 1x3 vector or 3 separate arguments')
end
phi   = x(1);
theta = x(2);
psi   = x(3);

y     = R1phi(phi)*R2theta(theta)*R3psi(psi);
end