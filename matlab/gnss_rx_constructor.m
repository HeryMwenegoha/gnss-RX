% Hery A Mwenegoha copyright (c) 2020

function gnssObj = gnss_rx_constructor(varargin)
% This a constructor function to allow running the GNSS measurement
% simulator in closed-loop mode. This is more representative of how the
% real-time implementation would work. 

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
config.odr_Hz      = 1;

% epoch format for 2 receivers mounted on the aircraft
fmt = struct('PRN',cell(1,32),'C1C',[],'L1C',[],'D1C',[],'S1C',[],'LLI',[]);
RXa = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);
RXb = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);

% Loop to extract configs from varargins
if nargin > 1
    for id = 1:nargin
        if ischar(varargin{id})
            switch(varargin{id})
                % number of runs
                case {'nruns', 'noruns', 'runs', 'run'}
                    if isnumeric(varargin{id+1})
                        config.no_runs = varargin{id+1};
                    end
                    
                case {'file', 'datafile', 'data'}
                    %user motion file
                    if ischar(varargin{id+1})
                        config.userFile = varargin{id+1};
                    end
                    
                case {'ephemfile', 'efile', 'ephemeris'}
                    %ephemeris file
                    if ischar(varargin{id+1})
                        config.ephemerisFile = varargin{id+1};
                    end
                    
                case {'odr', 'odr_Hz', 'odrHz', 'ODR'}
                    %output data rate
                    if isnumeric(varargin{id+1})
                        config.odr_Hz = varargin{id+1};
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

%% Initialise delays for our receiver
dtRx_s      = 1/config.odr_Hz;
gnssObj.Rxa = RXbias(dtRx_s);
gnssObj.Ir  = Iono(dtRx_s);
gnssObj.Tr  = Tropo(dtRx_s);
gnssObj.Tha = RXthermal(dtRx_s);
gnssObj.Thb = RXthermal(dtRx_s);
for index=1:32
    gnssObj.Mpa(index) = Multipath(dtRx_s);% multipath class
    gnssObj.Na(index)  = Nr;               % integer ambiguity class
end

%% Read IGS file
fprintf('-- reading Ephemeris -- \n');
[SV_IGS,IonoCoefficients] =readEphemeris(config.ephemerisFile);

%% Iono coefficients used here are perturbed by 10 percent
% user gets clean coefficients
fprintf('-- perturbing iono coefficients -- \n');
IonoCoeff.alpha  = IonoCoefficients.alpha + ...
    IonoCoefficients.alpha.*0.1.*randn(1,4);
IonoCoeff.beta   = IonoCoefficients.beta + ...
    IonoCoefficients.beta.*0.1.*randn(1,4);

% Store ionocofficients
gnssObj.IonoCoeff = IonoCoeff;

% Print Sample SpaceVehicle
SV_IGS(25);

% Number of Satellites loaded
nSats  = length(SV_IGS);

% Create a cell structure
SV  = struct('ID', cell(nSats,1));

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

% Store SV ephemeris parameters
gnssObj.svGps = SV;

% Set my experiment time as the max TOC time
fprintf('-- set exp start time from max TOC -- \n');
time_flght_start = SV(idxMax).TOC_g_time;

% Configure DOY
config.DOY = time_flght_start.g_doy;

config.start_g_time_s = time_flght_start;

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

% Store the errors to include
gnssObj.include = include;

% Store the configurations
gnssObj.config = config;
end