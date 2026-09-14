% Copyright © Hery A Mwenegoha copyright 2024 - 2026

classdef Gnss < handle
    properties
        mask_angle    = 2;     % mask angle [degrees]
        baseline1     = 1;     % baseline
        information   = [];    % information
        rcv2_enbaled  = false; % enable second receiver
        odr_Hz        = 1;     % output data rate
        incRnd        = false; % include random noise
        incIono       = true;  % include iono delays
        incTropo      = true;  % include tropo delays
        incMp         = true;  % include multipath
        incTh         = true;  % include thermal noise
        incRxClk      = true;  % include receiver clock errors e.g. if
                               % false then the clockBias and clockDrift
                               % will be zero
        incTxClk      = true;  % include satellite clock errors.
                               % is not using the mask
        incNcyc       = true;  % include integer number of cycles in
                               % carrier-phase measurements
        ephemerisFile =...     % sample ephemeris file
            fullfile('IGS','ABMF00GLP_R_20190611300_01H_GN.rnx');
        no_runs       = 1;     % [NOT USED] number of runs
        IonoCoeff     =...
            struct('alpha',[],...
                   'beta', []);
        gpsStartTime_s g_time = g_time(datetime(2019,03,02,13,00,00));
    end

    properties(SetAccess = private)
        RxaClk RXbias;         % Receiver clock model
        Ir     Iono;           % Ionospheric delay model
        Tr     Tropo;          % Tropospheric delat model
        Tha    RXthermal;      % Thermal noise model
        Thb    RXthermal;      % Thermal noise model
        Mpa    Multipath;      % Multipath model
        Na     Nr;             % Carrier-phase integer ambiguity
        svGps  struct;
        obs    struct;
    end

    methods
        function obj = Gnss(opts)
            % Constructor for the Gnss object
            arguments
                opts.?Gnss
            end

            for prop = string(fieldnames(opts))'
                obj.(prop) = opts.(prop);
            end

            % Some printouts
            fprintf('...Settings:\n');
            fprintf('...ephemerisFile\t:\t%1s\n', obj.ephemerisFile);

            % initialise classes
            dtRx_s    = 1/obj.odr_Hz;
            obj.RxaClk= RXbias(dtRx_s);
            obj.Ir    = Iono(dtRx_s);
            obj.Tr    = Tropo(dtRx_s);
            obj.Tha   = RXthermal(dtRx_s);
            obj.Thb   = RXthermal(dtRx_s);
            for index=1:32
                obj.Mpa(index) = Multipath(dtRx_s);
                obj.Na(index)  = Nr;
            end

            % Read IGS file
            disp('...reading Ephemeris');
            [EphStruct, IonoCoefficients] = readEphemeris(obj.ephemerisFile);

            % iono-coefficients for the klobuchar model
            if isempty(obj.IonoCoeff.alpha)
                fprintf('...perturbing iono coefficients \n');
                obj.IonoCoeff.alpha  = IonoCoefficients.alpha + ...
                    IonoCoefficients.alpha.*0.1.*randn(1,4);
                obj.IonoCoeff.beta   = IonoCoefficients.beta + ...
                    IonoCoefficients.beta.*0.1.*randn(1,4);
            end
            disp(obj.IonoCoeff);

            obj.svGps = Gnss.createSv(EphStruct);

            % Set my experiment time as the max TOC time
            if obj.gpsStartTime_s == g_time(datetime(2019,03,02,13,00,00))
                fprintf('...set exp start time from max TOC \n');
                tocValues = nan(1, length(obj.svGps));
                hasToc    = ~arrayfun(@(s) isempty(s.TOC), obj.svGps);
                tocValues(hasToc) = [obj.svGps(hasToc).TOC];
                [~,idxMax]   = max(tocValues);
                expStartTime = obj.svGps(idxMax).TOC_g_time;
                obj.gpsStartTime_s = expStartTime;
                disp(obj.gpsStartTime_s);
                fprintf('...tow_s: \t %.2f \n',obj.gpsStartTime_s.g_sow);
            end

            % Print DOY
            fprintf('...DOY: \t %u \n',obj.gpsStartTime_s.g_doy);
        end

        function [rawxSoln] = update(obj, current_SOW, posEcef, velEcef, rpy)
            % This function is called every navEpoch to compute the raw GNSS
            % observables.
            % Inputs:
            %   obj - the gnssObj structure containing information for the
            %   simulated receiver.
            %   current_SOW - the current time given in SOW.
            %   posEcef - current position (m) in ECEF coordinates - 3x1
            %   velEcef -  current velocity (m) in ECEF coordinate frame - 3x1
            %   Rpy - current attitude in Euler angles 3 x 1
            % Outputs:
            %   gnssObj - updated gnssObj. The different delays are updated through
            %   time
            %   rawxSoln - the output receiver solution containing raw GNSS observables
            %   returned for each succesful satellite
            arguments
                obj Gnss
                current_SOW (1,1) double
                posEcef (3,1) double
                velEcef (3,1) double
                rpy     (3,1) double
            end

            % Current epoch g_time
            t = current_SOW;

            % initialise rawxSoln
            rawxSoln = Gnss.emptyObs();

            % Get the current receiver position from the Navigation engine
            r_ea_e1   = posEcef; % 3x1
            v_ea_e1   = velEcef; % 3x1
            att.roll  = rpy(1);  % 1x1
            att.ptch  = rpy(2);  % 1x1
            att.yaw   = rpy(3);  % 1x1

            % TODO: we definitely want a way of checking this per navEpoch. But the
            % Nav Engine epoch (incoming pose solution) and gnss-RX epochs are
            % different for our case.
            % Number of Satellites in SV struct
            nSats  = length(obj.svGps);

            % Get the DOY
            sowLapsed = current_SOW - obj.gpsStartTime_s.g_sow;
            gpsTime   = obj.gpsStartTime_s + seconds(sowLapsed);
            DOY       = gpsTime.g_doy;

            % Pack the config expected by sat.PR - same for every
            % satellite this epoch, so build it once
            config.mask_angle = obj.mask_angle;
            config.DOY        = DOY;
            config.incRnd     = obj.incRnd;
            config.incIono    = obj.incIono;
            config.incTropo   = obj.incTropo;
            config.incMp      = obj.incMp;
            config.incTh      = obj.incTh;
            config.incRxClk   = obj.incRxClk;
            config.incTxClk   = obj.incTxClk;
            config.incNcyc    = obj.incNcyc;

            % calculate satellite positions at time of transmission - t
            for iPrn=1:nSats
                if isempty(obj.svGps(iPrn).ID)
                    continue;
                end

                prn = obj.svGps(iPrn).ID;

                % Get classes
                class.config    = config;
                class.IonoCoeff = obj.IonoCoeff;
                class.Ir        = obj.Ir;
                class.Tr        = obj.Tr;
                class.Mp        = obj.Mpa(prn);
                class.Th        = obj.Tha;
                class.Rx        = obj.RxaClk;
                class.Nn        = obj.Na(prn);
                class.r_e       = r_ea_e1;
                class.v_e       = v_ea_e1;
                class.att       = att;

                % Compute measurements
                [gT_r,C1C,L1C,D1C,LLI,cn0, satPos, satVel, dTs, ddTs]=sat.PR(obj.svGps, prn,t,class);

                % Store measurements if OK
                if ~isnan(C1C)
                    % Receiver data holder
                    rawxSoln.gTr          = gT_r;
                    rawxSoln.SOW          = gT_r*1;       % can be presented ok
                    rawxSoln.DOY          = DOY;

                    % receiver clockBias and drift terms
                    rawxSoln.bias_m      = obj.RxaClk.delay*wgs84.c;
                    rawxSoln.drift_mps   = obj.RxaClk.ddelay*wgs84.c;

                    % Note: slightly redundant information about the
                    %     : driving solution (truth solution this epoch)
                    rawxSoln.pos_ecef    = r_ea_e1;
                    rawxSoln.vel_ecef    = v_ea_e1;

                    % note: svg is indexed by prn and it is pre-allocated
                    % to the expected PRN size per constellation
                    % svg - 36, sve - 36
                    rawxSoln.svg(prn).PRN = prn;
                    rawxSoln.svg(prn).C1C = C1C;
                    rawxSoln.svg(prn).L1C = L1C;
                    rawxSoln.svg(prn).D1C = D1C;
                    rawxSoln.svg(prn).S1C = [];
                    rawxSoln.svg(prn).cn0 = cn0;
                    rawxSoln.svg(prn).LLI = LLI;
                    rawxSoln.svg(prn).pos = satPos;
                    rawxSoln.svg(prn).vel = satVel;
                    rawxSoln.svg(prn).clkBias_m = dTs*wgs84.c_light;
                    rawxSoln.svg(prn).clkDrift_mps = ddTs*wgs84.c_light;
                end
            end

            % persistent obs struct
            obj.obs = rawxSoln;

            % RX-A CLOCK MODEL
            obj.RxaClk.common();

            % IONO-GM Process
            obj.Ir.common();

            % TROPO-GM Process
            obj.Tr.common();
        end
    end

    methods(Static)
        function obs = emptyObs()
            % Used to initialise empty obs struct
            fmtGps = struct('PRN',cell(1,32),'C1C',[],'L1C',[],'D1C',[],'S1C',[],'LLI',[], 'cn0', []);
            fmtGal = struct('PRN',cell(1,36),'C1C',[],'L1C',[],'D1C',[],'S1C',[],'LLI',[], 'cn0', []);
            obs    = struct('SOW',[],'DOY',[], 'gTr', [], 'svg',fmtGps,'sve',fmtGal);
        end

        function SV = createSv(EphStruct)
            % Function that takes an ephemeris structure creates a space
            % vehicle structure.
            arguments
                EphStruct struct
            end

            % Number of Satellites loaded
            nSats = length(EphStruct);

            % Create a cell structure
            SV= struct('ID', cell(nSats,1));

            % Populate SV structure
            for i = 1:nSats
                SV(i).ID         = EphStruct(i).ID;
                SV(i).Health     = EphStruct(i).HEALTH;
                SV(i).eo         = EphStruct(i).e0;       %
                SV(i).toe        = EphStruct(i).TOE;      % ofGPSWEEK
                SV(i).Io         = EphStruct(i).IO;       % Orbital Inclination
                SV(i).OMEGAd_dot = EphStruct(i).OMEGADOT; % RateofRightAscencion [rad/sec]
                SV(i).a          =(EphStruct(i).SQRT_A).^2;
                SV(i).OMEGAo     = EphStruct(i).OMEGA0;   % RightAscenatWeek
                SV(i).w          = EphStruct(i).omega;    % ArgumentofPerigee
                SV(i).Mo         = EphStruct(i).M0;       % MeanAnom;    % Mean anomally at t=0      -- [change] at perigee
                SV(i).dn         = EphStruct(i).DELTA_N;  % Mean motion correction
                SV(i).Crs        = EphStruct(i).CRS;      % Sine   harmonic  radius   correction term --[m]
                SV(i).Crc        = EphStruct(i).CRC;      % Cosine harmonic  radius   correction term --[m]
                SV(i).Cus        = EphStruct(i).CUS;      % Sine   harmonic  argofLat correction      --[rad]
                SV(i).Cuc        = EphStruct(i).CUC;      % Cosine harmonic  argofLat correction      --[rad]
                SV(i).Cis        = EphStruct(i).CIS;      % Sine   harmonic  inclina  correction      --[rad]
                SV(i).Cic        = EphStruct(i).CIC;      % Cosine harmonic  inclina  correction      --[rad]
                %SV(i).satgroup  = mySat;
                SV(i).Id_dot     = EphStruct(i).IDOT;     % Not in Yuna [rad/sec]
                SV(i).TOC        = EphStruct(i).TOC;      % Reference Epoch of Clock data in secOfGPSWeek
                SV(i).A0         = EphStruct(i).A0;       % Sv clock bias [seconds]
                SV(i).A1         = EphStruct(i).A1;       % Sv clock drift [seconds/seconds]
                SV(i).A2         = EphStruct(i).A2;       % Sv clock drift rate [seconds/seconds2]
                SV(i).TGD        = EphStruct(i).TGD;      % Group delay [s]
                SV(i).TOC_g_time = EphStruct(i).TOC_g_time;      % Group delay [s]
            end

            % Sanity Checks
            toe = [SV.toe];
            toeDiffMinMax = max(toe) - min(toe);

            toc = [SV.TOC];
            tocDiffMinMax = max(toc) - min(toc);

            if toeDiffMinMax > 86400
                warning("... minMaxToe differ by %u seconds",toeDiffMinMax);
            end

            if tocDiffMinMax > 86400
                warning("...minMaxToc differ by %u seconds",tocDiffMinMax);
            end

            if any((toc - toe) > 0)
                warning("...Some Tocs differ from Toe");
            end
        end
    end
end