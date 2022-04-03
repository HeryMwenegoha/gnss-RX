% Hery A Mwenegoha copyright (c) 2020

function varargout = ecef2geodetic(xe, ye, ze, refEllipsoid)
    if strcmpi(refEllipsoid, 'wgs84')
        [varargout{1},varargout{2}, varargout{3}]=wgs84.ecef2geod_rad(xe, ye, ze);
    end
end