% Hery A Mwenegoha copyright (c) 2020

function [xNorth,yEast,zDown] = ecef2ned(X,Y,Z,lat0,lon0,h0,spheroid)
% Get local ECEF-coordinates
[Xe0, Ye0, Ze0]=geodetic2ecef(lat0, lon0, h0, spheroid);

% Delta ECEFs
dX    = X - Xe0;
dY    = Y - Ye0;
dZ    = Z - Ze0;

% ECEF-to-NED
R     = Ren(lat0, lon0);
dNED  = R(:,1)*dX(:).' + R(:,2)*dY(:).' + R(:,3)*dZ(:).';

if size(X,1) == 1
    xNorth=dNED(1,:);
    yEast =dNED(2,:);
    zDown =dNED(3,:);
else
    xNorth=dNED(1,:).';
    yEast =dNED(2,:).';
    zDown =dNED(3,:).';    
end
end