% Author    : Hery A Mwenegoha copyright (c) 2018
% version
%   current : 1.01    19/05/2019
%   first   : 1.00    01/07/2018
%
% Description
%   WGS84 Model
%
% syntax
%   y=wgs84
%
% input
%
%
% output
%   
%
% Contact
%   email   : hery_amani@yahoo.com

classdef wgs84 < handle
    properties (Constant)
        radius_equator     = 6378137;
        radius_poles       = 6356752.31424518;   % [m]
        earth_eccentricity = 0.0818191908426215; % [units]
        mu                 = 3.9860050e14;       % [m^3/s^2] earts gravitational constant - 3.986004418e14;   
        omega_ie           = 7.2921151467e-5;    % [rad/s]   earth-rotation rate
        hgt_gps_sat        = 20180000;           % [m] average height of gps sats
        TILT               = deg2rad(23*0); % Used in animation only
        c                  = 299792458;          % [m/s] speedOfLight in freespace
        Flattening         = 1 / 298.257223563;  % Earth flattening e^2 = 2f-f^2
        
        % For backward compatibility only
        semimajor          = 6378137; 
        eccentricity       = 0.0818191908426215;
        c_light            = 299792458;
    end
    
    methods
        % gravity  model
        
        % magnetic model
    end
    
    methods(Static)
        function y=e
            % Eccentricity
            y=wgs84.eccentricity;
        end
        
        function y=e2
            % Eccentricty-squared
            y=wgs84.f*(2-wgs84.f);
        end
        
        function y = Re
            % Equatorial radius
            y=wgs84.radius_equator;
        end
        
        function y = Rpo
            % Polar radius
            y=wgs84.radius_poles;
        end
        
        function y = A
            % Semimajor=Equatorial axis
            y=wgs84.semimajor;
        end
        
        function y = f
            % Flattening
            y=wgs84.Flattening;
        end
        
        function RM=RM(lat)
            % Meridian radius of curvature
            RM = (wgs84.A.*(1 - wgs84.e.^2))./(1 - wgs84.e.^2.*(sin(lat)).^2).^(3/2); 
        end
        
        function RP=RP(lat)
            % Prime Vertical radius of curvature
            % Transverse radius of curvature
            RP = wgs84.A./sqrt(1 - wgs84.e^2.*(sin(lat)).^2);                       
        end
        
        function RE=RE(lat)
            % Transverse radius of curvature: E-W motion
            RE=wgs84.RP(lat);
        end
         
        function varargout = geod2ecef_rad(lat, lon, h, varargin)
            % Geodetic to ECEF frame
            e2 = wgs84.f*(2-wgs84.f);
            v  = wgs84.RE(lat);
            x  = (v+h).*cos(lat).*cos(lon);
            y  = (v+h).*cos(lat).*sin(lon);
            z  = (v.*(1-e2)+h).*sin(lat);
            varargout{1} = x;
            varargout{2} = y;
            varargout{3} = z;
        end
        
        function varargout = ecef2geod_rad(xe, ye, ze)
            % J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates 
            % to geodetic coordinates," IEEE Transactions on Aerospace and 
            % Electronic Systems, vol. 30, pp. 957-961, 1994.
            % Heikkinen's Formula
            lambda =  atan2(ye, xe);
            
            x   = xe;
            y   = ye;
            z   = ze;

            r   = sqrt(x.^2+y.^2);
            a   = wgs84.A;
            e   = wgs84.e;
            esq = e*e;
            e2  = esq;
            b   = a*sqrt(1-esq);

            % Heikkinen's Formula
            if norm([x,y,x]) > 52000
                Esq = (a^2-b^2)/b^2;
                F   = 54*b^2*z.^2;
                G   = r.^2 + (1-esq).*z.^2 - esq*(a^2-b^2);
                c   = e^4.*F.*r.^2./G.^3;
                s   = (1 + c + (c.^2+2*c).^(1/2)).^(1/3);
                P   = F./(3.*(s+1./s+1).^2.*G.^2);
                Q   = (1+2*e^4.*P).^(1/2);
                r0  = -(P.*e^2.*r)./(1+Q)+(a^2./2*(1+1./Q)-P.*(1-e^2).*z.^2./(Q.*(1+Q))-P.*r.^2./2).^(1/2);
                U   = ((r-e^2.*r0).^2+z.^2).^(1/2);
                V   = ((r-e^2.*r0).^2+(1-e^2).*z.^2).^(1/2);
                z0  = b^2.*z./(a.*V);
                hb  = U.*(1 - b^2./(a.*V));
                mu  = atan((z+Esq.*z0)./r);
            else
                % Direct bowrings formula if below the threshold
                p   = sqrt(x.*x+y.*y);
                r   = sqrt(p.*p+z.*z);
                u   = atan(b.*z.*(1+e.*b./r)./(a.*p));
                mu  = atan((z+e.*b.*sin(u).^3)./(p-e2.*a.*cos(u).^3));
                v   = a./sqrt(1-e2.*sin(mu).^2);
                hb  = p.*cos(mu)+z.*sin(mu)-a*a./v;
            end
            
            varargout{1} = mu;
            varargout{2} = lambda;
            varargout{3} = hb;
        end   
    end
end