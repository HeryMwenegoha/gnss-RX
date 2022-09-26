% Author    : Hery A Mwenegoha (C) UoN 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   RX Oscillator drift and bias - RX Clock Model adopted from
%   Robert,Introduction to Random Signals.
% Functions
%   delay   -   common receiver bias in seconds for each signal
%   common  -   GM process loop
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

classdef RXbias < handle
    properties(Constant)
         delta_t = 1;      % [1 seconds sampling interval]
         ho      = 2e-20;  % 2e-20 gives us 8.9876e-04 m2/s   Tiny noise on the  2e-21
         h_2     = 2e-20;  % 2e-20 gives us 0.0355  m2/s3, where 0.04m2/s3 [groves] original 3e-24
    end
    
    properties(Access=private)
        bias;
        drift;
        offset;
        offset_rate;
    end
    
    methods 
        function this = RXbias
            this.offset      = 300000;  %3e5 [m]
            this.offset_rate = 20;      %20  [m/s]
            this.bias        = randn*this.offset/wgs84.c;
            this.drift       = randn*this.offset_rate/wgs84.c;%sqrt(2*math.PI^2*this.h_2*this.delta_t);
        end
        
        function y=delay(this)
           y=this.bias; % [seconds]
        end
        
        function y=ddelay(this)
            y=this.drift;
        end
        
        function common(this, input)
            % spectral densities for frequency and phase
            Sf = this.ho/2;
            Sg = 2*PI^2*this.h_2;
            
            qf  = Sf * this.delta_t;  % frequency
            qg  = Sg * this.delta_t;  % phase
            
            this.bias  = this.bias  +...
                this.drift*this.delta_t +...
                randn*sqrt(qf);
            this.drift = this.drift + randn*sqrt(qg);
        end
    end
    
end