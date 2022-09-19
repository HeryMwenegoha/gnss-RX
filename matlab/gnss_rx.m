% Hery A Mwenegoha copyright (c) 2020
% Lite version of the Satellite File
% Raw observables output with receiver timestamp

function DataOut = gnss_rx(varargin)
% Reset the random number seed for repeatability
rng(1);

% Error models to include
include.Random = false;
include.Iono   = true;
include.Tropo  = true;
include.MP     = true;
include.Th     = true;
include.Rx     = true;
include.Tx     = true;
include.N      = true;

% User Configs | Also through varargins
config.mask_angle  = 15;
config.baseline1   = 1;
config.information = [];%information;
config.rcv2_enbaled=false;
config.no_runs     = 1;
config.userFile    = ['simulator' filesep 'simple2d-03-04.mat'];
%'Trajectory\SIMDATA_AEROSONDE_100Hz_[f_ib_omega_ib]_TECS2_SND2_04_04_19_wen.mat';
config.ephemerisFile=['IGS' filesep 'ABMF00GLP_R_20190611300_01H_GN.rnx'];

% epoch format for 2 receivers mounted on the aircraft
fmt = struct('PRN',cell(1,32),'C1C',[],'L1C',[],'D1C',[],'S1C',[],'LLI',[]);
RXa = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);
RXb = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);

% Loop to extract configs from varargins
if nargin > 1
    for id = 1:nargin
        if ischar(varargin{id})
            switch(varargin{id})
                case {'nruns', 'noruns', 'runs', 'run'}
                    if isnumeric(varargin{id+1})
                        config.no_runs = varargin{id+1};
                    end
                %user motion file    
                case {'file', 'datafile', 'data'}
                    if ischar(varargin{id+1})
                        config.userFile = varargin{id+1};
                    end
                %ephemeris file
                case {'ephemfile', 'efile', 'ephemeris'}
                    if ischar(varargin{id+1})
                        config.ephemerisFile = varargin{id+1};
                    end                  
            end
        end
    end
elseif nargin == 1
    config.no_runs = varargin{1};
end

% Some printouts
fprintf('Settings:\n');
fprintf('Runs\t\t:\t%1u\n',config.no_runs);
fprintf('userMotionFile\t:\t%1s\n',config.userFile);
fprintf('ephemerisFile\t:\t%1s\n',config.ephemerisFile);

TimeLapse   = 0; % Total Time Lapse monte
%% Initialise delays
%for idx=1:config.no_runs

Rxa = RXbias;
Rxb = RXbias;
Ir  = Iono;
Tr  = Tropo;
Tha = RXthermal;
Thb = RXthermal;
for index=1:32
    Mpa(index)=Multipath;    % multipath class
    Mpb(index)=Multipath;    % multipath class
    Na(index) =Nr;           % integer ambiguity class
    Nb(index) =Nr;           % integer ambiguity class
end

% Read IGS file
[SV_IGS,IonoCoefficients] =readEphemeris(config.ephemerisFile);

% Iono coefficients used here are perturbed by 10 percent
% user gets clean coefficients
IonoCoeff.alpha  = IonoCoefficients.alpha + ...
    IonoCoefficients.alpha.*0.1.*randn(1,4);
IonoCoeff.beta   = IonoCoefficients.beta + ...
    IonoCoefficients.beta.*0.1.*randn(1,4);

% Print Sample SpaceVehicle
SV_IGS(25);

% Number of Satellites loaded
nSats  = length(SV_IGS);

% Create a cell structure
SV     = struct('ID', cell(nSats,1));

minTOC = 604800;
maxTOC = 0;
idxMax = 0;
idxMin = 0;

% Populate SV structure
for i = 1:nSats
    SV(i).ID         = SV_IGS(i).ID;
    SV(i).Health     = SV_IGS(i).HEALTH;
    SV(i).eo         = SV_IGS(i).e0;       %
    SV(i).toe        = SV_IGS(i).TOE;      % ofGPSWEEK
    SV(i).Io         = SV_IGS(i).IO;       % Orbital Inclination
    SV(i).OMEGAd_dot = SV_IGS(i).OMEGADOT; % RateofRightAscencion [rad/sec]
    SV(i).a          =(SV_IGS(i).SQRT_A).^2;
    SV(i).OMEGAo     = SV_IGS(i).OMEGA0;   % RightAscenatWeek
    SV(i).w          = SV_IGS(i).omega;    % ArgumentofPerigee
    SV(i).Mo         = SV_IGS(i).M0;       % MeanAnom;    % Mean anomally at t=0      -- [change] at perigee
    SV(i).dn         = SV_IGS(i).DELTA_N;  % Mean motion correction
    SV(i).Crs        = SV_IGS(i).CRS;      % Sine   harmonic  radius   correction term --[m]
    SV(i).Crc        = SV_IGS(i).CRC;      % Cosine harmonic  radius   correction term --[m]
    SV(i).Cus        = SV_IGS(i).CUS;      % Sine   harmonic  argofLat correction      --[rad]
    SV(i).Cuc        = SV_IGS(i).CUC;      % Cosine harmonic  argofLat correction      --[rad]
    SV(i).Cis        = SV_IGS(i).CIS;      % Sine   harmonic  inclina  correction      --[rad]
    SV(i).Cic        = SV_IGS(i).CIC;      % Cosine harmonic  inclina  correction      --[rad]
    %SV(i).satgroup  = mySat;
    SV(i).Id_dot     = SV_IGS(i).IDOT;     % Not in Yuna [rad/sec]
    SV(i).TOC        = SV_IGS(i).TOC;      % Reference Epoch of Clock data in secOfGPSWeek
    SV(i).A0         = SV_IGS(i).A0;       % Sv clock bias [seconds]
    SV(i).A1         = SV_IGS(i).A1;       % Sv clock drift [seconds/seconds]
    SV(i).A2         = SV_IGS(i).A2;       % Sv clock drift rate [seconds/seconds2]
    SV(i).TGD        = SV_IGS(i).TGD;      % Group delay [s]
    SV(i).TOC_g_time = SV_IGS(i).TOC_g_time;      % Group delay [s]
    
    if ~isempty(SV(i).TOC_g_time)
        if minTOC > SV(i).TOC
            minTOC = SV(i).TOC;
            idxMin = i;
        end
        
        if maxTOC < SV(i).TOC
            maxTOC = SV(i).TOC;
            idxMax = i;
        end
    end
end

% Set my experiment time as the max TOC time
time_flght_start = SV(idxMax).TOC_g_time;

% Configure DOY
config.DOY = time_flght_start.g_doy;

% Sanity Checks
TOE = [SV.toe];
TOE_Diff = max(TOE) - min(TOE);

TOC = [SV.TOC];
TOC_Diff = max(TOC) - min(TOC);

if TOE_Diff > 86400
    fprintf("---------------------------\n");
    warning("The TOE differ by %u seconds",TOE_Diff);
    fprintf("---------------------------\n");
end

if TOC_Diff > 86400
    fprintf("---------------------------\n");
    warning("The TOC differ by %u seconds",TOC_Diff);
    fprintf("---------------------------\n");
end

if any((TOC - TOE) > 0)
    fprintf("---------------------------\n");
    warning("Some TOCs differ from TOEs");
    fprintf("---------------------------\n");    
end

%% Load trajectory from user motion file
[mprofile, maxSeconds] = loadTrajectory(config.userFile);
epoch = 0;

% Store Ionospheric delay
Iono_delay    = zeros(0,0);
Tropo_delay   = zeros(0,0);

% DT_step is driven by the input user motion. 
% TODO: use navEpochPeriod to define this option externally
DT_step = round(mean(diff(mprofile.time)));
if DT_step >= 30
    maxSeconds = 3600 * 24;
end

% what time are we flying : GPS time of WEEK
simulationTime_sec = 0; % Total current lapse
last_perc   = 0;
for t=time_flght_start.g_sow:DT_step:time_flght_start.g_sow+maxSeconds-1 %3600*24
    tic;
    
    % Current epoch
    epoch    = epoch + 1;
    if epoch == 1
        if DT_step > 1
            fprintf('epoch and DT_step not in-sync  \n');
        end
    end
    
    Percentage = floor((epoch)/maxSeconds*100);
    if mod(Percentage,10) == 0 && Percentage ~= last_perc
        fprintf('Monte-completion%10.2f%s\n',Percentage, '%');
    end
    last_perc  = Percentage;
    
    % Get the DOY
    sow_lapsed= t-time_flght_start.g_sow;
    g_t2      = time_flght_start + seconds(sow_lapsed);
    DOY       = g_t2.g_doy;
    
    % current true receiver positions
    r_ea_e1   = mprofile.r_ea_e(:,epoch);
    v_ea_e1   = mprofile.v_ea_e(:,epoch);
    att.roll  = mprofile.roll(epoch);
    att.ptch  = mprofile.ptch(epoch);
    att.yaw   = mprofile.yaw(epoch);
    
    % calculate satellite positions at time of transmission - t
    increment = 0;
    for im=1:nSats
        %id
        if isempty(SV(im).ID) == true
            continue;
        end
        
        id = SV(im).ID;
         
        % Pseudorange
        class.config    = config;
        class.include   = include;
        class.IonoCoeff = IonoCoeff;
        class.Ir        = Ir;
        class.Tr        = Tr;
        class.Mp        = Mpa(id);
        class.Th        = Tha;
        class.Rx        = Rxa;
        class.Nn        = Na(id);
        class.r_e       = r_ea_e1;
        class.v_e       = v_ea_e1;
        class.att       = att;
        
        [gT_r,C1C,L1C,D1C,LLI]=sat.PR(SV, id,t,class);
        
        if ~isnan(C1C)
            % Receiver data holder
            % -- LLI - flag
            RXa(epoch).gTr         = gT_r;
            RXa(epoch).SOW         = gT_r*1; % can be presented ok
            RXa(epoch).DOY         = DOY;
            RXa(epoch).svg(id).PRN = SV(id).ID;
            RXa(epoch).svg(id).C1C = C1C;
            RXa(epoch).svg(id).L1C = L1C;
            RXa(epoch).svg(id).D1C = D1C;
            RXa(epoch).svg(id).S1C = [];
            RXa(epoch).svg(id).LLI = LLI;
            RXa(epoch).bias_m      = Rxa.delay*wgs84.c;
            RXa(epoch).drift_mps   = Rxa.ddelay*wgs84.c;
            RXa(epoch).pos_ecef    = r_ea_e1;
            RXa(epoch).vel_ecef    = v_ea_e1;
            
            increment              = increment + 1;
            %t_record_true(epoch)        = t;
            %t_record(epoch)             = gT_r*1;
            %dt_record(epoch)            = gT_r.f;
            %t_clk_bias(epoch)           = Rxa.delay;
            %t_clk_drift(epoch)          = Rxa.ddelay;
            %Iono_delay(epoch,increment) = delay.C1C.Ir;
            %Tropo_delay(epoch,increment)= delay.C1C.Tr;
        end
        
    end
    %}
    % RX-A CLOCK MODEL
    Rxa.common;
    
    % RX-B CLOCK MODEL
    Rxb.common;
    
    % IONO-GM Process
    Ir.common;
    
    % TROPO-GM Process
    Tr.common;
    
    % Zenith Iono delay
    Iono_delay(epoch)=Ir.residual;
    
    % Housekeeping
    Interval   = toc;
    simulationTime_sec    = simulationTime_sec   + Interval;  % Total current lapse
    TimeLapse  = TimeLapse + Interval;  % Total Time Lapse monte    
end
fprintf('Monte-Epoch: %1u , SimTime: %10.2f s , TimeLapse: %10.2f s \n', epoch, simulationTime_sec, TimeLapse);

RXa(epoch+1:end)=[];
RXb(epoch+1:end)=[];

%% Record data per simulation run
DataOut.Iono_delay      = Iono_delay; % Only the residual
DataOut.IonoCoefficients = IonoCoefficients;
DataOut.Function  = 'gnss_rx';
DataOut.include   = include;
DataOut.RXa       = RXa;
DataOut.RXb       = RXb;
DataOut.SVg_e     = SV; % GPS ephemeris only | TODO: svG_e, svE_e, svR_e
end