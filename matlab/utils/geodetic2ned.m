% Hery A Mwenegoha copyright (c) 2020

function [xNorth,yEast,zDown] = geodetic2ned(lat,lon,h,lat0,lon0,h0,spheroid)
% Get ECEF coordinates
[X, Y, Z]=geodetic2ecef(lat, lon, h, spheroid);

% Get Origin ECEF-coordinates
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