% Author    : Hery A Mwenegoha (C) 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   Ionospheric Modeling(common) and Return Delay per Satellite(delay)[s]
% Functions
%   common -   common ionospheric delay
%   delay  -  signal specific delay [s]
% syntax
%   common(this)
%     y=delay(this, input, params)
%   input
%     EL_degrees  : Elevation         [deg]
%     Az_degrees  : Azimuth           [deg]
%     lat_degrees : user latitude     [deg]
%     lon_degrees : user longitude    [deg]
%     tGPSweek    : timeinGPSWeek     [secs]
%   params
%     alpha0 = params.alpha(1);
%     alpha1 = params.alpha(2);
%     alpha2 = params.alpha(3);
%     alpha3 = params.alpha(4);
%
%     beta0  = params.beta(1);
%     beta1  = params.beta(2);
%     beta2  = params.beta(3);
%     beta3  = params.beta(4);
%   output
%      y   : time delay/advance due to ionospheric effect
classdef Iono < handle
    % IONOSPHERIC MODEL
    properties(Constant)
        % average sigma values
        % 0.45m    Ephemeris Error
        % 1.00m    Satellite Clock
        % 0.20m    Troposphere Error
        % 3-4m     Troposphere Error []
        
        sigma     = 1.7; % 0.7 best zenith value
        Tau       = 1800;
    end
    
    properties(Access = private)
        % Internal use only
        qk;
        residual_;
        residual_old;
    end
    
    properties
        delta_t;
    end
    
    methods
        % Constructor
        function this = Iono(delta_t)
            arguments
                delta_t double = 1.0; % sample rate [s]
            end
            this.delta_t   = delta_t;
            this.qk        = this.sigma.^2*(1 - exp(-2*this.delta_t/this.Tau));
            this.residual_ = this.sigma*randn;
            this.residual_old = 0;
        end
        
        % Common ionospheric residual in metres
        function common(this)
            this.residual_old = this.residual_;
            this.residual_  = ...
                (this.residual_.*exp(-this.delta_t/Iono.Tau) +...
                randn*sqrt(this.qk));
        end
        
        % Ionospheric residual common [m]
        function y=residual(this)
            y = this.residual_;
        end
        
        % Rate of change of ionospheric delay [s]
        function y=ddelay(this)
            y = ((this.residual_ - this.residual_old)./wgs84.c)./this.delta_t;
        end
        
        % Total delay in [seconds]
        function y=delay(this, input, params)
            [ts,F] = Iono.kb_seconds(input, params);
            y = ts + F*this.residual_/wgs84.c;
        end
    end % METHODS - PUBLIC
    
    methods(Static)
        % Estimate Ionospheric time delay [s]
        function y=estimate(input, params)
            [y,F] = Iono.kb_seconds(input, params);
        end
        
        % klobuchar model : kb_s : klobuchar_seconds
        % Estimate ionospheric time delay Y [s] and slant factor F
        function [y,F] = kb_seconds(input, params)
            if nargin >= 1
                if isa(input, 'struct')
                    EL_degrees  = input.EL_degrees;
                    Az_degrees  = input.Az_degrees;
                    lat_degrees = input.lat_degrees;
                    lon_degrees = input.lon_degrees;
                    tGPSweek    = input.tGPSweek;
                end
            end
            
            % important function
            deg2semi = 1.0/180;
            semi2rad = PI;
            
            if nargin < 2
                error('Iono::delya -- less parameters input ');
            else
                alpha0 = params.alpha(1);
                alpha1 = params.alpha(2);
                alpha2 = params.alpha(3);
                alpha3 = params.alpha(4);
                
                beta0  = params.beta(1);
                beta1  = params.beta(2);
                beta2  = params.beta(3);
                beta3  = params.beta(4);
            end
            
            % Earth centered angle from user to IPP
            % [semicircles]
            Psi_pp = 0.0137./(EL_degrees.*deg2semi + 0.11) - 0.022;
            
            % latitude of ionospheric pierce point (IPP) from user latitude
            % [semicirclces]
            lat_ipp = lat_degrees*deg2semi +...
                Psi_pp.*cos(deg2rad(Az_degrees));
            
            if lat_ipp > 0.416
                lat_ipp = 0.416;
            end
            if lat_ipp < -0.416
                lat_ipp = -0.416;
            end
            
            % longitude of ionospheric pierce point from user's longitude
            % [semicircles]
            lon_ipp = lon_degrees*deg2semi + ...
                Psi_pp.*sin(deg2rad(Az_degrees))./cos(lat_ipp*semi2rad);
            
            % Geomagnetic latitude
            % [semicircles]
            lat_geom = lat_ipp + 0.064.*cos((lon_ipp - 1.617)*semi2rad);
            
            % local time from IPP longitude and actual GPS time [s]
            % local time of the day at IPP
            t   = 43200.*lon_ipp + tGPSweek;
            t   = mod(t, 86400);
            
            if t > 86400
                t = t - 86400;
            end
            if t < 0
                t = t + 86400;
            end
            
            % Slant factor
            F = 1 + 16.*(0.53 - EL_degrees*deg2semi).^3;
            
            % Amplitude of ionospheric delay [seconds]
            AMP = alpha0              + ...
                alpha1.*lat_geom      +...
                alpha2.*lat_geom.^2   +...
                alpha3.*lat_geom.^3;
            if AMP < 0
                AMP = 0;
            end
            
            % Period of delay [seconds]
            PER = beta0               + ...
                beta1.*lat_geom       +...
                beta2.*lat_geom.^2    +...
                beta3.*lat_geom.^3;
            if PER < 72000
                PER = 72000;
            end
            
            % Phase of ionospheric delay
            XI = 2.*PI.*(t-50400)./PER; % [radians]
            
            % Ionospheric time delay [seconds].
            if abs(XI) < 1.57
                Ir_s=F.*(5.*10^-9 + AMP.*(1-XI.^2/2+XI.^4/24));
            else
                Ir_s=F.*(5.*10^-9);
            end
            y=Ir_s;
        end
        
        % klobuchar model : slant-factor
        function F = kb_sf(EL_degrees)
            % important function
            deg2semi = 1.0/180;
            
            % Slant factor
            F = 1 + 16.*(0.53 - EL_degrees*deg2semi).^3;
        end
    end % METHODS - STATIC
end