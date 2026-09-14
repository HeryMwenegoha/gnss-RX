% Author    : Hery A Mwenegoha (C) 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   GM process Multipath Modeling per satellite signal. The time const
%   ant is randomly initialised per signal
% Functions
%   delay         -     Return signal delay s]
% syntax
%   y=delay(this, input, params)
% input
%   EL_degrees  : Elevation         [deg]
%   Az_degrees  : Azimuth           [deg]
%   lat_degrees : user latitude     [deg]
%   lon_degrees : user longitude    [deg]
%   tGPSweek    : timeinGPSWeek     [secs]
% output
%   y   : time delay/advance due to ionospheric effect
classdef Multipath < handle
    properties(Constant)
        C0        = 0.47;  %[m]
        C1        = 0.78;  %[m]
        C2        = 20.92; %[deg]
        
        % low elevation satellites shows 2mm of sigma with this
        C0_ca     = 0.0002;%0.0015;%[m]
        C1_ca     = 0.001;%0.006; %[m]
        C2_ca     = 35;    %[deg]
    end
    
    properties(Access = private)
        residual;
        residual_ca;
        residual_ca_old;
    end
    
    properties
        Tau;
        delta_t;
    end
    
    methods
        % Constructor
        function this = Multipath(delta_t)
            arguments
                delta_t (1,1) double  =  1; % sample rate [s]
            end
            this.Tau          = randi([3,40],1);
            this.residual     = 0;
            this.residual_ca  = 0;
            this.residual_ca_old = 0;
            this.delta_t = delta_t;
        end
        
        % Pseudorange multipath [seconds]
        function Mp = delay(this, input, params)
            El              = input.EL_degrees;
            % Standard deviation of driving noise
            sigma_mp        = this.C0 + this.C1.*exp(-El./this.C2);
            qk              = sigma_mp.^2*(1 - exp(-2*this.delta_t/this.Tau));
            this.residual   = this.residual.*exp(-this.delta_t/this.Tau) + randn*sqrt(qk); %[m]
            Mp              = this.residual./wgs84.c; %[s]
        end
        
        % Carrier phase mutlipath [s]
        function Mp = delay_ca(this, input, params)
            El      = input.EL_degrees;
            
            % Standard deviation of driving noise
            sigma_mp = this.C0_ca + this.C1_ca.*exp(-El./this.C2_ca);
            qk       = sigma_mp.^2*(1 - exp(-2*this.delta_t/this.Tau));
            this.residual_ca_old=this.residual_ca;
            this.residual_ca= this.residual_ca.*exp(-this.delta_t/this.Tau)...
                + randn*sqrt(qk);                 %[m]
            Mp       = this.residual_ca./wgs84.c; %[s]
        end
        
        % No units - Meant for the doppler measurements
        function y = ddelay(this)
            y = ((this.residual_ca - this.residual_ca_old)./wgs84.c)./this.delta_t;
        end
    end
end