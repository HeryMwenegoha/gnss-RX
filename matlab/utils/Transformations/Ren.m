% Hery A Mwenegoha copyright (c) 2020

function y=Ren(x)
% ECEF-to-NED transformation
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