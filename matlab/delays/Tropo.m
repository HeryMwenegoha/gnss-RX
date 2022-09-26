% Author    : Hery A Mwenegoha (C) 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   Tropospheric Modeling(common) and Return Delay per Satellite(delay)[s]
% Functions
%   common             -    common tropospheric delay
%   delay              -    signal specific delay [s]
%   EGNOS              -    EGNOS interpolation tables
%   M_black_and_eisner -    black and eisner mapping functions
% syntax
%   common(this)
%   y=delay(this, input, params)
% input:
%   EL_degrees  : Elevation         [deg]
%   Az_degrees  : Azimuth           [deg]
%   lat_degrees : user latitude     [deg]
%   lon_degrees : user longitude    [deg]
%   tGPSweek    : timeinGPSWeek     [secs]
% output
%   y   : time delay/advance due to ionospheric effect
classdef Tropo < handle
    properties(Constant)
        sigma     = 0.1; %0.20m    Troposphere Error
        delta_t   =  1;
        Tau       = 1800;
    end
    
    properties(Access = private)
        % Internal use only
        qk;
        residual_;
        residual_old;
    end
    
    methods
        function this=Tropo
            this.qk        = this.sigma.^2*(1 - exp(-2*this.delta_t/this.Tau));
            this.residual_ = this.sigma*randn;
            this.residual_old = 0;
        end
        
        function y = residual(this)
            y =  this.residual_;
        end
        
        function common(this)
            % common residual process to all satellites at every epoch
            this.residual_old = this.residual_;
            this.residual_=(this.residual_.*exp(-this.delta_t/this.Tau)...
                + randn*sqrt(this.qk));
        end
        
        function y = ddelay(this)
            y = ((this.residual_ - this.residual_old)./wgs84.c)./this.delta_t;
        end
        
        function Tr = delay(this, input, params)
            % Tropo delay [s]
            % satellite specific delay in seconds
            if isa(input, 'struct')
                El          = input.EL_degrees;  % user elevation  [deg]
                lat_degrees = input.lat_degrees; % user latitude   [deg]
                DOY         = input.DOY;         % Day of the Year [day]
                H           = input.Hb; % user altitude above mean sea level - orthometric [m]
            else
                errordlg('input not struct');
                return;
            end
            
            egnos = Tropo.EGNOS(lat_degrees, DOY);
            
            % constants
            g    = 9.80665;              %[m/s2]
            
            T    = egnos.temperature;    % Sea level Temperature [k]
            
            Beta = egnos.Beta;           % Temperature lapse rate [K/m]
            
            Rd   = 287.054;              % Gas Constant [J/Kg/K]
            
            nu   = egnos.nu;             % water vapour lapse rate
            
            k1   = 77.604;               % K/mbar
            
            k2   = 382000;               % K2/mbar
            
            P    = egnos.pressure;       % mbar
            
            gm   = 9.784;                % m/s^2
            
            %https://en.wikipedia.org/wiki/Vapour_pressure_of_water
            e    = egnos.water_vapour;  % mbars at mean sea level - water vapour pressure
            
            ZTDd0 = (10^-6).*k1.*Rd.*P/gm;
            
            ZTDw0 = (10^-6).*k2.*Rd./(gm*(nu+1) - Beta.*Rd).*e./T;
            
            % Zenith Hydrostatic Delay [m]
            ZTDd = ZTDd0.*(1 - Beta.*H./T).^(g./(Rd.*Beta));
            
            % Zenith Wet Delay [m]
            ZTDw = ZTDw0.*(1 - Beta.*H./T).^((nu+1).*g./(Rd.*Beta) - 1);
            
            % Zenith Total Delay [m]
            ZTDt = ZTDd + ZTDw + this.residual;
            
            % Slant Tropospheric Delay [sec]
            Tr   = ZTDt.*Tropo.M_black_and_eisner(El)'./wgs84.c;
        end
    end
    
    methods(Static)
        function Tr = estimate(input, params)
            % satellite specific delay in seconds using EGNOS model
            if isa(input, 'struct')
                El          = input.EL_degrees;  % user elevation  [deg]
                lat_degrees = input.lat_degrees; % user latitude   [deg]
                DOY         = input.DOY;         % Day of the Year [day]
                H           = input.Hb;          % user altitude above mean sea level - orthometric [m]
            else
                errordlg('input not struct');
                return;
            end
            
            egnos = Tropo.EGNOS(lat_degrees, DOY);
            
            % constants
            g    = 9.80665;              %[m/s2]
            
            T    = egnos.temperature;    % Sea level Temperature [k]
            
            Beta = egnos.Beta;           % Temperature lapse rate [K/m]
            
            Rd   = 287.054;              % Gas Constant [J/Kg/K]
            
            nu   = egnos.nu;             % water vapour lapse rate
            
            k1   = 77.604;               % K/mbar
            
            k2   = 382000;               % K2/mbar
            
            P    = egnos.pressure;       % mbar
            
            gm   = 9.784;                % m/s^2
            
            %https://en.wikipedia.org/wiki/Vapour_pressure_of_water
            e    = egnos.water_vapour;  % mbars at mean sea level - water vapour pressure
            
            ZTDd0 = (10^-6).*k1.*Rd.*P/gm;
            
            ZTDw0 = (10^-6).*k2.*Rd./(gm*(nu+1) - Beta.*Rd).*e./T;
            
            % Zenith Hydrostatic Delay [m]
            ZTDd = ZTDd0.*(1 - Beta.*H./T).^(g./(Rd.*Beta));
            
            % Zenith Wet Delay [m]
            ZTDw = ZTDw0.*(1 - Beta.*H./T).^((nu+1).*g./(Rd.*Beta) - 1);
            
            % Zenith Total Delay [m]
            ZTDt = ZTDd + ZTDw;
            
            % Slant Tropospheric Delay [sec]
            Tr   = ZTDt.*Tropo.M_black_and_eisner(El)'./wgs84.c;
        end
        
        function Tr = estimate_Saastamoinen(input, percRelHumidity)
            % Simple Tropospheric model
            % satellite specific delay in seconds
            if isa(input, 'struct')
                El          = input.EL_degrees;  % user elevation  [deg]
                lat_degrees = input.lat_degrees; % user latitude   [deg]
                DOY         = input.DOY;         % Day of the Year [day]
                H           = input.Hb; % user altitude above mean sea level - orthometric [m]
            else
                errordlg('input not struct');
                return;
            end
            
            h    = H;  % altitude amsl
            if nargin >= 2
                hrel = percRelHumidity;
            else
                hrel = 70; % relative humidity at 70% approximate
            end
            P    = 1013.15*(1 - 2.2557 * 10^-5 * h)^5.2568;
            T    = 15.0 - 6.5*10^-3*h + 273.15;
            e    = 6.108 * exp((17.15*T - 4684.0)/(T-38.45))* hrel./100;
            
            Z    = pi/2 - deg2rad(El); % Zenith angle
            Tr   = (0.002277./cos(Z).*(P + (1255/T + 0.05)*e - tan(Z).^2)./wgs84.c).';
        end
        
        function y=EGNOS(lat_degrees, DOY)
            % EGNOS Tropospheric Model
            % lat_degrees : latitude in degrees of User
            % DOY         : Day of the Year of the Flight
            Latitude_Chart    = [15      30      45       60       75];
            
            % Average Values
            Pressure_Chart    = [1013.25 1017.25 1015.75  1011.75  1013.00]; % mbar
            Temperature_Chart = [299.65  294.15  283.15   272.15   263.65];  % K
            water_vapor_Chart = [26.31   21.79   11.66    6.78     4.11];    % mbar
            Beta_chart        = [6.30    6.05    5.58     5.39     4.53];    % mK/m
            nu_chart          = [2.77    3.15    2.57     1.81     1.55];    % []
            
            % Seasonal Variation
            dPressure_Chart    = [0	-3.75	-2.25	-1.75	-0.50]; % mbar
            dTemperature_Chart = [0	 7.00	11.00	15.00	14.50]; % K
            dwater_vapor_Chart = [0	 8.85	 7.24	 5.36	 3.39]; % mbar
            dBeta_chart        = [0	 0.25	 0.32	 0.81	 0.62]; % mK/m
            dnu_chart          = [0	 0.33	 0.46	 0.74	 0.30]; % []
            
            if lat_degrees >= 0
                DOYmin = 28;
            else
                DOYmin = 211;
            end
            
            lat = abs(lat_degrees);
            if lat < Latitude_Chart(1)
                lat = Latitude_Chart(1);
            elseif lat > Latitude_Chart(end)
                lat = Latitude_Chart(end);
            end
            
            pressure_mean     = interp1(Latitude_Chart, Pressure_Chart,     lat, 'linear');
            temperature_mean  = interp1(Latitude_Chart, Temperature_Chart,  lat, 'linear');
            water_vapour_mean = interp1(Latitude_Chart, water_vapor_Chart,  lat, 'linear');
            Beta_mean    = interp1(Latitude_Chart, Beta_chart,         lat, 'linear');
            nu_mean      = interp1(Latitude_Chart, nu_chart,           lat, 'linear');
            
            dpressure    = interp1(Latitude_Chart, dPressure_Chart,     lat, 'linear');
            dtemperature = interp1(Latitude_Chart, dTemperature_Chart,  lat, 'linear');
            dwater_vapour= interp1(Latitude_Chart, dwater_vapor_Chart,  lat, 'linear');
            dBeta        = interp1(Latitude_Chart, dBeta_chart,         lat, 'linear');
            dnu          = interp1(Latitude_Chart, dnu_chart,           lat, 'linear');
            
            pressure       = pressure_mean - dpressure.*cos(2.*pi.*(DOY - DOYmin)./365.25);
            temperature    = temperature_mean - dtemperature.*cos(2.*pi.*(DOY - DOYmin)./365.25);
            water_vapour   = water_vapour_mean - dwater_vapour.*cos(2.*pi.*(DOY - DOYmin)./365.25);
            Beta           = Beta_mean - dBeta.*cos(2.*pi.*(DOY - DOYmin)./365.25);
            nu             = nu_mean - dnu.*cos(2.*pi.*(DOY - DOYmin)./365.25);
            
            y.pressure     = pressure;
            y.temperature  = temperature;
            y.water_vapour = water_vapour;
            y.Beta         = Beta*10^-3;
            y.nu           = nu;
        end
        
        function mapp=M_black_and_eisner(El_degrees)
            % black and eisner mapping function
            % El_degrees : Elevation in degrees
            El   = El_degrees;
            gamma= 1 - 0.015*max(0,4 - El).^2;
            mapp = 1.001./((0.002001 + (sin(deg2rad(El))).^2).^(1/2)).*gamma;
        end
    end
end