% Author    : Hery A Mwenegoha (C) UoN 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   Receiver Thermal Noise as a function of CNO
% Functions
%   delay   -   Return Delay per observable(delay)[s]
%   CNO_    -   Carrier-to-Noise density ratio
%   sigma   -   SD of code or carrier based on CNO
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
classdef RXthermal < handle
    properties(Constant)
        sigma_ca = 0.0007;
        delta_t  = 1;
    end
    
    properties (Access = public)
        % access the returned SD in [m] for code and carrier
        sigma_;
    end
    
    properties(Access = private)
        L1;
        antx_gain_pattern;
        res_ca;
        res_ca_old;
    end
    
    methods
        function this= RXthermal
            % Constructor
            this.L1.C0   = 0.05; % [m]
            this.L1.C1   = 1.05; % [m]
            this.L1.C2   = 28.0; %[dB-Hz]
            this.L1.C3   = 8.00; %[dB-Hz]
            
            El_in_degrees    = [ 0  15 30 65  90];
            Antenna_Gain_dBi = [-6  -1  3  9  12];
            this.antx_gain_pattern.El = El_in_degrees;
            this.antx_gain_pattern.Ga = Antenna_Gain_dBi;
            this.res_ca_old = 0;
            this.res_ca  = 0;
        end
        
        function y=delay(this,input)
            % Satellite specific delay in seconds
            El_in_deg   = input.EL_degrees;
            if El_in_deg     > max(this.antx_gain_pattern.El)
                El_in_deg = max(this.antx_gain_pattern.El);
            elseif El_in_deg < min(this.antx_gain_pattern.El)
                El_in_deg    = min(this.antx_gain_pattern.El);
            end
            
            cno         = CNO_(this,this.antx_gain_pattern, El_in_deg);
            sigma_all   = sigma(this,cno);          % [m]
            this.sigma_ = sigma_all;
            y           = sigma_all.code * randn * 1/wgs84.c;  % [sec]
        end
        
        % carrier phase delay [s]
        function y=delay_ca(this,input)
            El_in_deg   = input.EL_degrees;
            if El_in_deg     > max(this.antx_gain_pattern.El)
                El_in_deg = max(this.antx_gain_pattern.El);
            elseif El_in_deg < min(this.antx_gain_pattern.El)
                El_in_deg    = min(this.antx_gain_pattern.El);
            end
            cno      = CNO_(this,this.antx_gain_pattern, El_in_deg);
            sigma_all = sigma(this,cno);          % [m]
            
            this.res_ca_old = this.res_ca;
            this.res_ca = sigma_all.carrier * randn * 1/wgs84.c; % [sec]
            y           = this.res_ca;  % [sec]
        end
        
        % No units - Meant for the doppler measurements
        function y = ddelay(this)
            y = (this.res_ca - this.res_ca_old)./this.delta_t;
        end
        
        function y=CNO_(this,gain_pattern, El_in_deg)
            % Signal to Noise Ratio
            El_in_degrees    = gain_pattern.El;
            Antenna_Gain_dBi = gain_pattern.Ga;
            Pr  = -158.5;  %[dBW]Power of GNSS signal at receiver position
            Ga  =  interp1(El_in_degrees,Antenna_Gain_dBi, El_in_deg);   %[dB] Receiving antenna gain
            No  = -201;    %[dB-Hz] Noise density ratio
            SD  =    1;    %[dB]
            y   = Pr + Ga - No  + SD*randn; %
        end
        
        function y=sigma(this,CNO_dB_Hz)
            % Calculate the standard deviations of the code and carrier
            CNO           = CNO_dB_Hz;     % carrier-to-noise density ratio between receiver and satellite [dB-Hz]
            sigma_e       = this.L1.C0 + this.L1.C1.*exp(-(CNO - this.L1.C2)./this.L1.C3);
            sigma_carrier = this.sigma_ca;
            y.code        = sigma_e;       % [m]
            y.carrier     = sigma_carrier; % [-]
        end
    end
    
    methods(Static)
        function plt()
            fig  = figure('Position', [201 114 610 501]);
            El_in_degrees    = [ 0  15 30 65  90];
            Antenna_Gain_dBi = [-6  -1  3  9  10];
            El_in_deg        = 0:90;
            Pr  = -158.5;  %[dBW]Power of GNSS signal at receiver position
            Ga  =  interp1(El_in_degrees,Antenna_Gain_dBi, El_in_deg);   %[dB] Receiving antenna gain
            No  = -201;    %[dB-Hz] Noise density ratio
            SD  =    1;    %[dB]
            %y   = Pr + Ga - No  + SD*randn; %
            y   = Ga;
            plot(El_in_deg, y);
            
            xlabel('Elevation [deg]');
            ylabel('Antenna Gain [dBi]');
            grid on;
            box on;
            set(gca, 'FontSize', 14, 'XLim', [0,90],'TickLength',[0,0]);
        end
    end
end
