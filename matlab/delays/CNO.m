% By Hery A Mwenegoha (C) 2020

classdef CNO < handle
    properties(Constant)
        El_in_degrees    = [ 0  15 30 65  90];
        Antenna_Gain_dBi = [-6  -1  3  9  12];
    end
    
    methods(Static)
        function c = cno(El_in_deg)
            % limits
            ip = El_in_deg > 90;
            if any(ip)
                El_in_deg(ip) = 90;
            end
            
            ip = El_in_deg < -90;
            if any(ip)
                El_in_deg(ip) = -90;
            end
            
            El_in_degrees    = CNO.El_in_degrees;
            Antenna_Gain_dBi = CNO.Antenna_Gain_dBi;
            
            % Power of GNSS signal at receiver position
            Pr  = -158.5;  %[dBW]
            
            % Receiver antenna gain interpolated to elevation
            Ga  =  interp1(El_in_degrees,...
                Antenna_Gain_dBi,...
                El_in_deg,'linear', 'extrap');%[dB]
            
            % Noise-Power-Density
            No  = -201;    %[dB-Hz]
            
            % Standard deviation of received noise
            SD  =    1;    %[dB]
            
            % output
            c   = Pr + Ga - No  + SD*randn;
        end
        
        function plot(El_in_deg)
            c = CNO.cno(El_in_deg);
            plot(El_in_deg, c);
            xlabel('Elevation [degrees]');
            ylabel('CNO [dB-Hz]');
            grid on;
            title('Carrier-to-Noise density');
        end
    end
end