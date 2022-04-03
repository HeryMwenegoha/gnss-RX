% Hery A Mwenegoha copyright (c) 2020

function y=Ren(varargin)
% ECEF-to-NED transformation
if nargin == 1 || nargin == 2
    x = [varargin{:}];
else
    error('input is either 1x2 vector or 2 separate arguments')
end
latRad = x(1);
lonRad = x(2);
n = [-sin(latRad)*cos(lonRad);...
     -sin(latRad)*sin(lonRad);...
      cos(latRad)];

e = [-sin(lonRad);
      cos(lonRad);
      0];

d = [-cos(latRad)*cos(lonRad);   %d = cross(n,e);
     -cos(latRad)*sin(lonRad);
     -sin(latRad)];

y = [n.'; e.'; d.'];
end