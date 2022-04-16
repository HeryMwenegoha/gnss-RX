% Author    : Hery A Mwenegoha (C) 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   GM process - Signal specific Interfrequency Bias with exponentially damped
%   standard deviation based on time [s]
% syntax
%   y=delay(this, input, params)
% input
%   EL_degrees  : Elevation         [deg]
%   Az_degrees  : Azimuth           [deg]
%   lat_degrees : user latitude     [deg]
%   lon_degrees : user longitude    [deg]
%   tGPSweek    : timeinGPSWeek     [secs]
% params
%   time [secs]
% output
%   y   : IFbias delay [s]
classdef IFbias < handle
    properties(Constant)
        sigma_to = 0.005;  % [m]
        C0       = 86400;  % [s]
        Tau      = 1;      % [s]
        delta_t  = 1;      % [s]
    end
    
    properties(Access=private)
        residual;
    end
    
    methods 
        function this= IFbias
            this.residual = 0;
        end
        
        function y=delay(this, input, params)
            time        = params.time;
            sigma_IF    = this.sigma_to * exp(-time/this.C0);
            
            qk            = sigma_IF.^2*(1 - exp(-2*this.delta_t/this.Tau));                       % [m^2]
            this.residual = this.residual.*exp(-this.delta_t/this.Tau) + randn*sqrt(qk); % [m]
            
            y=this.residual./wgs84.c;
        end
    end
end