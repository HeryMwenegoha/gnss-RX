% Hery A Mwenegoha copyright (c) 2020

function varargout = geodetic2ecef(lat, lon, h, refEllipsoid)
    if strcmpi(refEllipsoid, 'wgs84')
        [varargout{1}, varargout{2}, varargout{3}]=wgs84.geod2ecef_rad(lat,lon,h);
    end
end