% Hery A Mwenegoha copyright (c) 2020

classdef test_GnssRx < matlab.unittest.TestCase
    properties(Constant)
        plottingEnabled = false;
    end
    
    methods (Test)
        function test_Iono(self)
            Input.EL_degrees  = 20;
            Input.Az_degrees  = 0;
            Input.lat_degrees = 52;
            Input.lon_degrees = 0;
            Input.tGPSweek    = 568800;
            
            IonoCoeff.alpha = [9.2739e-09 0 -6.0652e-08 0];
            IonoCoeff.beta  = [9.7864e+04 0 -2.0642e+05 0];
            
            nMonte   = 50;
            epochMax = 100;
            resIonoBuff_m = nan(nMonte,epochMax);
            for iReal = 1:nMonte
                Ir = Iono();
                for k=1:epochMax
                    Ir.common();
                    resIono_m = Ir.residual();
                    resIonoBuff_m(iReal,k) = resIono_m;
                    %Ir_s = Ir.delay(Input, IonoCoeff);
                end
            end
            
            rmsVal = sqrt(sum(resIonoBuff_m.^2,1)./nMonte);
            
            if self.plottingEnabled
                figure('Name', 'Residual-Iono-Delay');
                hold on;
                plot(resIonoBuff_m.', 'Linewidth', 2);
                plot( rmsVal, '-k', 'LineWidth',2);
                plot(-rmsVal, '-k', 'LineWidth',2);
                hold off;
            end
            
            fprintf('...qualifications for %u monte carlo runs \n', nMonte);
            fprintf('...max epoch utilised      \t %u \n', epochMax);
            self.verifyLessThan(rmsVal, 3.0);
        end
    end
end