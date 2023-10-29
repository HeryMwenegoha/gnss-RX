% By Hery A Mwenegoha (C) 2020
% All satellite related utility
%   Calculate satellite position
%   Calculate satellite elevation
%   Calculate eccentric anomaly

classdef sat < handle
    methods(Static)
        function y=lambda1
            % L1 wavelength
            y = wgs84.c/sat.f1;
        end    
        
        function y=f1
            % L1 Frequency
            y=10.23*154e6;
        end      
        
        function y=f2
            % L2 frequency
            y=10.23*120e6;
        end
        
        function [r_es_e, varargout] = pos(SV, id, Tst, iter)
            % Satellite position
            if nargout == 1
                if nargin > 3
                   r_es_e  = sat.calculate_satpos(SV, id, Tst,iter); 
                else
                    r_es_e = sat.calculate_satpos(SV, id, Tst);
                end
            else
                if nargout == 2
                    if nargin > 3
                       [r_es_e, varargout{1}]= sat.calculate_satpos(SV, id, Tst,iter); 
                    else
                       [r_es_e, varargout{1}]= sat.calculate_satpos(SV, id, Tst);
                    end    
                elseif nargout == 3
                    if nargin > 3
                       [r_es_e, varargout{1}, varargout{2}]=...
                           sat.calculate_satpos(SV, id, Tst,iter); 
                    else
                       [r_es_e, varargout{1}, varargout{2}]=...
                           sat.calculate_satpos(SV, id, Tst);
                    end                    
                end
            end
        end
        
        function [r_es_e, varargout] = calculate_satpos(SVg_e, id, Tst, iter)
            % Satellite Position Calculation Given Time of Signal Transmission
            % SV  - Ephemeris Structure Object
            % id  - Index of Satellite in Ephemeris Object
            % Tst - Time of signal transmission [s] : estimated
            % iter - number of iterations for E-anom calculation
            iter_i = 1;
            
            if nargin > 3
                iter_i = iter;
            end
            
            % we need the ephemeris of this satellite
            SV_ = SVg_e(id);
            
            % Time objects
            if ~isa(Tst, 'gTime')
                gTst = gTime(Tst);
            else
                gTst = Tst;
            end
            gToe = gTime(SV_.toe);
            
            if SV_.a <= 0
                error('satpos error "A" <=0');
            end
            
            % calculated semi minor axis
            SV_.b = SV_.a * sqrt(1-SV_.eo^2);
            
            % toe : time of reference ephemeris
            %delta_t = Tst-SV_.toe;
            delta_t  = gTst-gToe;
            
            if (delta_t*1) > 3600*4
                %fprintf('...prn %u has outdated toe \n',id);
            end
            
            % week crossover check
            if abs(delta_t*1) > 302400
                % 604800 correction for week crossover
                if delta_t > 0
                    delta_t = delta_t - gTime(604800);
                elseif delta_t < 0
                    delta_t = delta_t + gTime(604800);
                end
                warning('Week crossover detected <calculate_satpos>');
            end
            
            % constants
            mu  = wgs84.mu;
            wis = sqrt(mu/SV_.a^3);
            % mean motion correction
            Delta_n = SV_.dn;
            % mean anomally
            M = SV_.Mo + (wis + Delta_n)*delta_t;
            
            % Eccentric anomally from keplers equation
            eo   = SV_.eo;
            E    =  M + (eo*sin(M))/(1-sin(M+eo)+sin(M));
            for ite=1:iter_i % 20 iterations cm accuracy
                E = M + eo*sin(E);
            end
            
            % True anomally
            v = atan2(sqrt(1-eo^2)*sin(E)/(1-eo*cos(E)),...
                (cos(E)-eo)/(1-eo*cos(E)));

            % Argument of latitude
            Theta = SV_.w+v; % Angle of satellite from right ascension
            
            % orbital radius and angle from nodal point
            % Crs - Sine   Harmonic corr. of orbital radius
            % Crc - Cosine Harmonic corr. of orbital radius
            a     = SV_.a;
            Crs   = SV_.Crs;
            Crc   = SV_.Crc;
            Cus   = SV_.Cus;
            Cuc   = SV_.Cuc;
                        
            % Orbital radius and argument of latitude
            r_os  = a*(1-eo*cos(E)) + Crs*sin(2*Theta) + Crc*cos(2*Theta);
            u_os  = Theta           + Cus*sin(2*Theta) + Cuc*cos(2*Theta);
            
            % orbital position
            x_os_o = r_os*cos(u_os);
            y_os_o = r_os*sin(u_os);
            z_os_o = 0;
            
            % Longitude of ascending node
            OMEGAo = SV_.OMEGAo;
            Io     = SV_.Io;
            Id_dot = SV_.Id_dot;
            Cis    = SV_.Cis; % 0; left this at zero for no apparent reason.
            Cic    = SV_.Cic; % 0; left this at zero for no apparent reason.
            TOE    = gTime(SV_.toe);
            OMEGAd_dot = SV_.OMEGAd_dot;
            OMEGA      = (OMEGAo - wgs84.omega_ie * (delta_t + TOE) +...
                OMEGAd_dot * delta_t);
            
            % wrapped to 360 degrees
            % OMEGA  = math.wrap_360(OMEGA*math.degreesf)*math.radiansf;
            Incl   = Io + Id_dot*delta_t + Cis*sin(2*Theta) + Cic*cos(2*Theta);
            
            E_dot     = (wis + Delta_n)/(1-eo*cos(E));
            Theta_dot = sin(v)/sin(E)*E_dot;
            r_os_dot  = (a*eo*sin(E))*E_dot +...
                2*(Crs*cos(2*Theta) - Crc*sin(2*Theta))*Theta_dot;
            u_os_dot  =(1+2*Cus*cos(2*Theta)-2*Cuc*sin(2*Theta))*Theta_dot;
            x_os_o_dot= r_os_dot*cos(u_os) - r_os*u_os_dot*sin(u_os);
            y_os_o_dot= r_os_dot*sin(u_os) + r_os*u_os_dot*cos(u_os);
            z_os_o_dot= 0;
            OMEGA_dot = OMEGAd_dot - wgs84.omega_ie;
            Incl_dot  = Id_dot+2*(Cis*cos(2*Theta)-Cic*sin(2*Theta))*Theta_dot;
            
            v_es_e=...
                [x_os_o_dot*cos(OMEGA)-y_os_o_dot*cos(Incl)*sin(OMEGA) + Incl_dot*y_os_o*sin(Incl)*sin(OMEGA)
                x_os_o_dot*sin(OMEGA)+y_os_o_dot*cos(Incl)*cos(OMEGA) - Incl_dot*y_os_o*sin(Incl)*cos(OMEGA)
                y_os_o_dot*sin(Incl)+Incl_dot*y_os_o*cos(Incl)]  + ...
                (wgs84.omega_ie - OMEGAd_dot)*...
                [x_os_o*sin(OMEGA) + y_os_o*cos(Incl)*cos(OMEGA)
                -x_os_o*cos(OMEGA) + y_os_o*cos(Incl)*sin(OMEGA)
                0];
            
            % Satellite Position Orbital frame
            r_os_o =[x_os_o;...
                y_os_o;...
                z_os_o];
            
            % Orbital frame to ECEF frame
            Roe    = [cos(OMEGA) -cos(Incl)*sin(OMEGA)  sin(Incl)*sin(OMEGA); ...
                      sin(OMEGA)  cos(Incl)*cos(OMEGA) -sin(Incl)*cos(OMEGA); ...
                      0           sin(Incl)             cos(Incl)];
            
            % ECEF frame position
            r_es_e = Roe * r_os_o;
            
            if nargout > 1
                varargout{1} = E;
                if nargout > 2
                    varargout{2} = v_es_e;
                end
            end
        end
        
        function y= elev(pos, r_es_e)
            % Satellite elevation
            y = sat.elevation(pos, r_es_e);
        end
        
        function [y,latRad, lonRad,h]= elev_ecef(r_ea_e, r_es_e)
            if size(r_es_e,2) > 1
                error('sat::elev_ecef::check r_es_e size');
            end
            
            if size(r_ea_e,2) > 1
                error('sat::elev_ecef::check r_ea_e size');              
            end
            
            % Geometric range in ECEF - sagnac correction
            GR    = r_es_e - r_ea_e ;
            GR_mag= sqrt(sum(GR.^2));%+sagna_correction;    
            
            % unit vector in ECEF coordinates
            GR_unit = GR./GR_mag;
            
            [latRad,lonRad,h]=wgs84.ecef2geod_rad(r_ea_e(1),r_ea_e(2),...
                                                  r_ea_e(3));
            
            n = [-sin(latRad)*cos(lonRad);...
                 -sin(latRad)*sin(lonRad);...
                  cos(latRad)];

            e = [-sin(lonRad);
                  cos(lonRad);
                  0];

            d = [-cos(latRad)*cos(lonRad);   %d = cross(n,e);
                 -cos(latRad)*sin(lonRad);
                 -sin(latRad)];
            u = -d;

            % Pack elevation and Azimuth
            y.EL  = asin(GR_unit.' * u);    
            y.Az  = atan2((GR_unit.'*e),(GR_unit.'*n)); % atan2 [-pi pi]
            y.LOS = GR_unit;
            y.lat = latRad;
            y.lon = lonRad;
            y.h   = h;
        end
        
        function y = elevation(pos, r_es_e)
            % uses user and satellite position
            % r_ea_e - 3x1 user position in ECEF
            % r_es_e - 3x1 satellite position
            [xe,ye,ze] = geodetic2ecef(wgs84Ellipsoid('metres'),...
                                       pos.lat,...
                                       pos.lon,...
                                       pos.hd,...
                                      'radians');    
            r_ea_e  = [xe; ye; ze];

            % sagnac correction effect for using ECEF frame
            sagna_correction = Sagnac.ecef(r_ea_e, r_es_e);

            % Geometric range in ECEF - sagnac correction
            GR    = r_es_e - r_ea_e ;
            GR_mag= sqrt(sum(GR.^2)) + sagna_correction;    

            % unit vector in ECEF coordinates
            GR_unit = GR./GR_mag;
            
            % user latitude and longitude needed
            latRad = pos.lat;
            lonRad = pos.lon;

            n = [-sin(latRad)*cos(lonRad);...
                 -sin(latRad)*sin(lonRad);...
                  cos(latRad)];

            e = [-sin(lonRad);
                  cos(lonRad);
                  0];

            d = [-cos(latRad)*cos(lonRad);   %d = cross(n,e);
                 -cos(latRad)*sin(lonRad);
                 -sin(latRad)];
            u = -d;

            % Pack elevation and Azimuth
            y.EL  = asin(GR_unit.' * u);    
            y.Az  = atan2((GR_unit.'*e),(GR_unit.'*n)); % atan2 [-pi pi]
            y.LOS = GR_unit;
        end       
        
        function Ek = E_anom(SV, id, Tst, loops)
            % Calculate Eccentric Anomaly
            if nargin > 3
               Ek = sat.E_anomaly(SV, id, Tst, loops); 
            else
               Ek = sat.E_anomaly(SV, id, Tst); 
            end
        end
        
        function Ek = E_anomaly(SV, id, Tst, loops)
            % Calculates an iterative Eccentric Anomaly based on:
            % SV  - Ephemeris Structure Object
            % id  - Satellite Index in Eph.Objec i.e. Not PRN
            % Tst - Satellite Time of Transmission <estimated>
            % loops - iterations for the Eccentric anomally calculation
            if nargin <= 3
                loops = 1;
            end
            
            % we need the ephemeris of this satellite
            SV_   = SV(id);
            
            if ~isa(Tst, 'gTime')
                gTst = gTime(Tst);
            else
                gTst = Tst;
            end
            gToe = gTime(SV_.toe);
            
            % calculated semi minor axis
            SV_.b = SV_.a * sqrt(1-SV_.eo^2);
            
            % our current time of transmission, i.e. transmit every second
            % difference between time of ephemeris and time of transmission
            %tst     = t;
            delta_t = gTst - gToe;
            
            % important
            % week crossover check
            if abs(delta_t*1) > 302400
                % 604800 correction for week crossover
                if delta_t > 0
                    delta_t = delta_t - gTime(604800);
                elseif delta_t < 0
                    delta_t = delta_t + gTime(604800);
                end
                fprintf('Warning : Week CrossOver Detected - check : calculate_Ek \n');
            end
            
            % Current Mean anomaly and Mean Motion Correction
            mu  = wgs84.mu;
            wis = sqrt(mu/SV_.a^3);
            M   = SV_.Mo + (wis + SV_.dn)*delta_t;
            
            % Eccentric anomally from keplers equation
            eo   = SV_.eo;
            E    =  M + (eo*sin(M))/(1-sin(M+eo)+sin(M));
            E2   = M;
            for ite=1:loops % 20 iterations cm accuracy
                E = M + eo*sin(E);
            end
            Ek = E;
        end
        
        function [gTr_o, C1C, L1C, D1C, LLI, cn0, r_es_e, v_es_e, dTs, ddTs]=PR(SV,id,trcv,class)
            % Simulate Pseudorange & Carrier Phase
            % Calculate actual transmission time
            % Return actual satellite position
            % Inputs:
            %   SV    - Satellite Ephemeris structure object
            %   id    - ID of the satellite to calculate measuremnent 1-32 for GPS
            %   trcv  - actual receiver time - noiseless
            %   class - Iono, Tropo,common receiver bias, integer amb, Navi
            %   [r-position ECEF,v-velocity ECEF,attitude], Iono-coefficients,m
            %   includes and other configs
            % Useable SV classes
            config    = class.config;
            include   = class.include;
            IonoCoeff = class.IonoCoeff;
            Ir        = class.Ir;
            Tr        = class.Tr;
            
            %if isfield(class, 'Mp')
            Mp        = class.Mp;
            %end
            
            Th        = class.Th;
            Rx        = class.Rx;
            
            %if isfield(class, 'Nn')
            Nn        = class.Nn;
            %end
            
            v_ea_e    = class.v_e;% 3x1
            r_ea_e    = class.r_e;
            att       = class.att;
            
            % rough intial estimate of sat position based on trcv
            [r_es_e0,Ekk] = sat.pos(SV, id, trcv, 10);
            r_es_e        = r_es_e0;
            
            % Day of the year
            DOY           = config.DOY;     
            GR_mag_last   = 0;
            dt_s_old      = 0;
            
            gTr           = gTime(trcv, 0);
            for i = 1:10
                % sagnac correction effect
                sagna_correction = Sagnac.ecef(r_ea_e, r_es_e);

                % Geometric range in ECEF - sagnac correction
                GR     = r_es_e - r_ea_e ;
                GR_mag = sqrt(sum(GR.^2)) + sagna_correction;                 
                
                % delta transmission
                % dt_s   = GR_mag/wgs84.c   + Ir_delay + Tr_delay;
                dt_s   = GR_mag/wgs84.c;
                
                % Tranmission time
                gTs    = gTr  - gTime(0, dt_s);

                delta_dT_s = dt_s - dt_s_old;
                dt_s_old   = dt_s;
                
                % Recalculate satellite position at new t
                [r_es_e, Ek, v_es_e] = sat.pos(SV, id, gTs, 23);

                if i ~= 1
                     delta_PR = GR_mag-GR_mag_last;
                     if abs(delta_PR) < 1e-6 % less than 1mm is acceptable
                         %Ts_obj = sat.gTime(trcv, -dt_s);
                         break;
                     elseif i == 10
                         fprintf('sat->actual:Ts did not converge \n');
                     end
                end
                GR_mag_last  = GR_mag;  % [GR_mag t_transmit]
            end
            
            % Geometric range 
            GR_unit= GR./GR_mag;  
            GR     = (gTr-gTs)*wgs84.c;
            
            % Input for Ir and Tr computations
            [y,latRad, lonRad,h]= sat.elev_ecef(r_ea_e, r_es_e);
            Ng                  = 0;       % Geoid Height [m]
            Hb                  = h - Ng;  % Othometric height [m]
            input.EL_degrees    = rad2deg(y.EL);
            input.Az_degrees    = rad2deg(y.Az);
            input.lat_degrees   = rad2deg(latRad);
            input.lon_degrees   = rad2deg(lonRad);
            input.Hb            = Hb;
            input.tGPSweek      = trcv;
            input.DOY           = DOY;               
            
            % Satellite clock delay
            [dTs,~,ddTs]=sat.dT_s(SV,id,gTs,1,'f1');
            %[dTs1,~,ddTs1]=sat.dT_s(SV,id,gTs+gTime(1));
            
            %{
            %if Nn.slp
            % Apparent Elevation
            e_us_n = [cos(y.EL)*cos(y.Az);
                      cos(y.EL)*sin(y.Az);
                     -sin(y.EL)];
            Rnb    = math.Rnb([att.roll, att.ptch, att.yaw]);
            zb     = Rnb(3,:);
            EL_apparent = -asin(zb * e_us_n);
            Nn.slip(EL_apparent*math.degreesf);
            %}
            
            % Calculated C1C delays
            Ir_delay    = Ir.delay(input,IonoCoeff)  * include.Iono;
            Tr_delay    = Tr.delay(input)            * include.Tropo; 
            Th_delay    = Th.delay(input)            * include.Th;  
            Mp_delay    = Mp.delay(input)            * include.MP;   
            Rx_delay    = Rx.delay                   * include.Rx;   
            Ts_delay    = dTs                        * include.Tx; 
            Rnd_delay   = (randn * 0.1)./wgs84.c     * include.Random;
            
            % Calculated carrier-phase delays
            Ir_delay_ca =  Ir.delay(input,IonoCoeff) * include.Iono;    
            Tr_delay_ca =  Tr.delay(input)           * include.Tropo;
            Th_delay_ca =  Th.delay_ca(input)        * include.Th;
            Mp_delay_ca =  Mp.delay_ca(input)        * include.MP;
            Rx_delay_ca =  Rx.delay                  * include.Rx;
            Ts_delay_ca =  dTs                       * include.Tx;
            Rnd_delay_ca=  randn * 0.0               * include.Random;
            N_integer   =  Nn.N                      * include.N; % cycles   
            
            % Delays for the Doppler Measurements
            Rx_delay_dd = Rx.ddelay                  * include.Rx; 
            Ir_delay_dd = Ir.ddelay                  * include.Iono;
            Tr_delay_dd = Tr.ddelay                  * include.Tropo;
            Mp_delay_dd = Mp.ddelay                  * include.MP; 
            Th_delay_dd = Th.ddelay                  * include.Th;
            
            % Corrupt Receiver Time | Common Error on Each PR Value
            gTr_o = gTr + gTime(0,Rx_delay);
            gTs_o = gTs + gTime(0,dTs);
            
            cn0 = CNO.cno(input.EL_degrees);
            if input.EL_degrees < config.mask_angle
                % cant produce measurements
                % fprintf('PRN %1u below mask angle \n', SV(id).ID);
                C1C = NaN;
                L1C = NaN;
                D1C = NaN;
                LLI = [];
                return;
            end

             sagnac_rate  = Sagnac.rate(r_ea_e, v_ea_e, r_es_e, v_es_e);
             
             %Pseudorange Measurement   [m]
             C1C = GR + (Rx_delay - Ts_delay).*wgs84.c + ...
                 (Ir_delay  +...
                 Tr_delay  +...
                 Mp_delay  +...
                 Th_delay  +...
                 Rnd_delay).*wgs84.c;
             
             % Carrier-phase Measurement [cycles]
             L1C = (GR + (Rx_delay_ca - Ts_delay_ca).*wgs84.c + ...
                 (-Ir_delay_ca  +...
                 Tr_delay_ca  +...
                 Mp_delay_ca*1+...
                 Th_delay_ca*1).*wgs84.c)./sat.lambda1 +...
                 N_integer*1;
             
             % Doppler frequency Measurement
             f1  =  wgs84.c/sat.lambda1;
             D1C = -f1./wgs84.c * ([v_es_e - v_ea_e].'*GR_unit +...
                 sagnac_rate                 +...
                 wgs84.c*Rx_delay_dd         -...
                 wgs84.c*ddTs                +...
                 wgs84.c*Ir_delay_dd         +...
                 wgs84.c*Tr_delay_dd         +...
                 wgs84.c*Th_delay_dd         +...
                 wgs84.c*Mp_delay_dd)         +...
                 randn*1.0; % 0.02
                              
            % Try to detect any cycle slips
            if Nn.detect == true
                LLI = 1;
                % fprintf('Cycle-Slip \n');
            else
                LLI = 0;
            end                                       
        end        
        
        function varargout=dT_s(SV,id,Ts, numIters_,tgd_correct_frequency)
            N        = nargin;
            dT_s_old = 0;
            T_s      = Ts;
            numIter  = 1;
            if nargin > 3
                numIter = numIters_;
            end
            b  =  0;
            fi =  sat.f1;
            if nargin > 4
                if ischar(tgd_correct_frequency)
                    switch(tgd_correct_frequency)
                        case {'F1', 'f1'}
                           fi =  sat.f1; 
                        case {'F2', 'f2'}
                           fi =  sat.f2; 
                        case {'F5', 'f5'}
                           fi = sat.f5;
                        case {'F3', 'f3'}
                           fi = sat.f3;
                        otherwise
                           error('HAM::sat.dT_s - frequency not defined');
                    end
                end
                b = sat.f1.^2/fi.^2;  
            end            
            
            for im=1:numIter
                Ek = sat.E_anom(SV,id,T_s,20);
                F  = -2*sqrt(wgs84.mu)/wgs84.c^2;
                dT_rel= F*...
                        SV(id).eo*...
                        sqrt(SV(id).a)*...
                        sin(Ek);
                dT_GD  = -b*SV(id).TGD;
                diff_Ts_Toc=T_s-gTime(SV(id).TOC);
                dT_s  = (SV(id).A0 +                ...
                         SV(id).A1.*diff_Ts_Toc    +...
                         SV(id).A2.*diff_Ts_Toc.^2 +...
                         dT_rel                    +...
                         dT_GD);
                delta_dT_s = dT_s - dT_s_old;
                dT_s_old   = dT_s;
                T_s        = T_s - gTime(0,delta_dT_s);
            end
            ddT_s  = SV(id).A1+2*SV(id).A2.*diff_Ts_Toc;            
                     
            if nargout > 1
                varargout{1} = dT_s;
                varargout{2} = T_s;
                if nargout > 2
                    varargout{3} = ddT_s;
                end
            else
                varargout{1} = dT_s;
            end
        end  
    end
end