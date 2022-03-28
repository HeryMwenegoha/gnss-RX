% Author: Hery A Mwenegoha (C) 2020 UoN-University of Nottingham NGI
% Writing this during the COVID-19 era 
% Lite version of the Satellite File
% Raw observables output with receiver time stamp
function DataOut = SatelliteLiteV3(varargin)
% Reset the random number seed
rng(1);

% Passed Information for the flight or day of event
% Should be able to decode this from the incoming msg
information.YearOfFlight    = 2019;
information.MonthOfFlight   = 03;
information.DayOfFlight     = 02;
information.GPSHourOfFlight = 14; 
information.GPSMinOfFlight  = 00;
information.GPSSecOfFlight  = 00;
time_flght_start            = g_time(...
                              datetime(information.YearOfFlight,...
                                       information.MonthOfFlight,...
                                       information.DayOfFlight,...
                                       information.GPSHourOfFlight,...
                                       information.GPSMinOfFlight,...
                                       information.GPSSecOfFlight),...
                              0);
information.GPSDayOfFlight = time_flght_start.g_doy;

% include
include.Random    = false;
include.Iono      = true;
include.Tropo     = true;
include.MP        = true;
include.Th        = true;
include.Rx        = true;
include.Tx        = true;
include.N         = true;

config.mask_angle  = 15;
config.baseline1   = 1;
config.information = information;
config.rcv2_enbaled=false;

% epoch format for 2 receivers mounted on the aircraft
fmt = struct('PRN',cell(1,32),'C1C',[],'L1C',[],'D1C',[],'S1C',[],'LLI',[]);
RXa = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);
RXb = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);

config.no_runs = 1;
config.file    = 'trajectorySimulator\2022-3-28-0-50-22-simple2d.mat';
                 %'Trajectory\SIMDATA_AEROSONDE_100Hz_[f_ib_omega_ib]_TECS2_SND2_04_04_19_wen.mat';
if nargin > 1
    for id = 1:nargin
        if ischar(varargin{id})
            switch(varargin{id})
                case {'nruns', 'noruns', 'runs', 'run'}
                    if isnumeric(varargin{id+1})
                        config.no_runs = varargin{id+1};
                    end
                    
                case {'file', 'datafile', 'data'}
                    if ischar(varargin{id+1})
                        config.file = varargin{id+1};
                    end
            end
        end
    end
elseif nargin == 1
    config.no_runs = varargin{1};
end

if isempty(config.file)
    error('specify <trajectory|file>');
end

fprintf('Settings:\n');
fprintf('Runs\t:\t%1u\n',config.no_runs);
fprintf('File\t:\t%1s\n',config.file);

return
TimeLapse   = 0; % Total Time Lapse monte
TimePerSVs  = 0;
for idx=1:config.no_runs
    % Initialise delays
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
        %Na(index).N
        %Nb(index).N
    end
    %return
    
    % Read IGS file
    [SV_IGS,IonoCoefficients] =...
    readEphemeris('IGS/ABMF00GLP_R_20190611300_01H_GN.rnx');
  

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
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Inputs From File 
    %TimeOfFlight = 568800;         % secondsofGPSWEEK - 1400hrs 
    maxSeconds   = loadTrajectory; % maxFlightTimeInSeconds
    epoch        = 0;
    %return

    % Storage 
    t_record_true = zeros(0,0);
    t_record      = zeros(0,0);
    dt_record     = zeros(0,0);
    t_clk_bias    = zeros(0,0);
    t_clk_drift   = zeros(0,0);
    Iono_delay    = zeros(0,0);
    %Iono_zenith_residual = zeros(0,0);
    Tropo_delay   = zeros(0,0);
    %Tropo_zenith_residual= zeros(0,0);
    %Ts_clk        = zeros(0,0);
    %GR_record     = zeros(0,0);
    %PR_record     = zeros(0,0);
    %L1C_record    = zeros(0,0);
    %Az_record     = zeros(0,0);
    %EL_record     = zeros(0,0);
    %ID_record     = zeros(0,0);

    %|-| update rate [seconds] - check if synchronized with trajectory %|-|200
    DT_step     = 1;
    if DT_step >= 30
       maxSeconds = 3600 * 24; 
    end
     
    for t=time_flght_start:DT_step:time_flght_start+seconds(maxSeconds)
        %t
    end
    
    % what time are we flying : GPS time of WEEK
    SimTime     = 0; % Total current lapse
    last_perc   = 0;
    for t=time_flght_start.g_sow:DT_step:time_flght_start.g_sow+maxSeconds %3600*24 
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
        r_ea_e1   = mprofile.r_ea_e1(epoch,:).';
        r_eb_e1   = mprofile.r_eb_e1(epoch,:).';
        v_ea_e1   = mprofile.v_ea_e(epoch,:).';
        v_eb_e1   = mprofile.v_eb_e(epoch,:).';
        att.roll  = mprofile.roll(epoch);
        att.ptch  = mprofile.ptch(epoch);
        att.yaw   = mprofile.yaw(epoch);
        
        % Same Time Frame
        %{
        ip  =~arrayfun(@(x)isempty(x.ID), SV); % non-empty SVs
        id  = [SV(ip).ID];                     % non-empty IDs
        class.config    = config;
        class.include   = include;
        class.IonoCoeff = IonoCoeff;
        class.Ir        = Ir;
        class.Tr        = Tr;
        %class.Mp        = Mpa(id);
        class.Th        = Tha;
        class.Rx        = Rxa;
        %class.Nn        = Na(id);
        class.r_e       = r_ea_e1;
        class.v_e       = v_ea_e1;
        class.att       = att;
        fun  = @(x_id) sat.PR(SV, x_id, t, class);
        [gT_r,C1C,L1C,D1C,LLI]=arrayfun(fun, id, 'UniformOutput', false);
        %}
        

        %
        % calculate satellite positions at time of transmission - t
        increment = 0;
        for im=1:nSats
            %id
            if isempty(SV(im).ID) == true
                continue;
            end
            
            id = SV(im).ID;

            % Receiver time using the gTime object
            %gTr       = gTime(t,0);
            %gTr_b     = gTime(t,0);
            
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
            
            %{
            class.config    = config;
            class.include   = include;
            class.IonoCoeff = IonoCoeff;
            class.Ir        = Ir;
            class.Tr        = Tr;
            class.Mp        = Mpb(id);
            class.Th        = Thb;
            class.Rx        = Rxb;
            class.Nn        = Nb(id);
            class.r_e       = r_eb_e1;
            class.v_e       = v_eb_e1;
            [gT_r_b,C1C_b,L1C_b,D1C_b,LLI_b]=sat.PR(SV, id,t,class);            
            
          
            if config.rcv2_enbaled == true  
                if ~isnan(C1C_b)
                    RXb(epoch).gTr         = gT_r_b;
                    RXb(epoch).SOW         = gT_r_b*1; % can bepresented ok
                    RXb(epoch).DOY         = DOY;
                    RXb(epoch).svg(id).PRN = SV(id).ID;
                    RXb(epoch).svg(id).C1C = C1C_b;     
                    RXb(epoch).svg(id).L1C = L1C_b;    % cycles
                    RXb(epoch).svg(id).D1C = D1C_b;
                    RXb(epoch).svg(id).S1C = [];
                    RXb(epoch).svg(id).LLI = LLI_b;   
                    RXb(epoch).bias_m      = Rxb.delay*wgs84.c;
                    RXb(epoch).drift_mps   = Rxb.ddelay*wgs84.c;
                    RXb(epoch).pos_ecef    = r_eb_e1;
                    RXb(epoch).vel_ecef    = v_eb_e1;
                end
            end
            %}
              
            
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
                
                increment                   = increment + 1;
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
        
        % Induce Cycle-Slip
        %if randi([1,32]) == 20
        %for index=1:32
        %Na(index) = Nr;% integer ambiguity class
        %Nb(index) = Nr;% integer ambiguity class
        %end
        %fprintf('Cycle - Slip \n');
        %end
        %fprintf('Epoch: %1u , Time: %1f s \n', epoch, t);
        
        Interval   = toc;
        SimTime    = SimTime   + Interval;  % Total current lapse
        TimeLapse  = TimeLapse + Interval;  % Total Time Lapse monte
        TimePerSVs = TimePerSVs + Interval;
        
        if nSats
            %fprintf('TimePer-32-SVs\t:\t%f\n',TimePerSVs);
            TimePerSVs = 0;
        end
    end
    fprintf('Monte-Epoch: %1u , SimTime: %10.2f s , TimeLapse: %10.2f s \n', idx, SimTime, TimeLapse);

    RXa(epoch+1:end)=[];
    RXb(epoch+1:end)=[];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record all the data per simulation run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DataOut(idx).t_record_true   = t_record_true;
    %DataOut(idx).t_record        = t_record;
    %DataOut(idx).dt_record       = dt_record;
    %DataOut(idx).t_clk_bias      = t_clk_bias;
    %DataOut(idx).t_clk_drift     = t_clk_drift;
    %DataOut(idx).Iono_residual   = Iono_zenith_residual; 
    DataOut(idx).Iono_delay       = Iono_delay;
    %DataOut(idx).Tropo_residual  = Tropo_zenith_residual;
    %DataOut(idx).Ts_clk          = Ts_clk;
    %DataOut(idx).Tropo_delay     = Tropo_delay; % getting negative values
    %DataOut(idx).PR_record       = PR_record;
    %DataOut(idx).L1C             = L1C_record;
    %DataOut(idx).Az_record       = Az_record;
    %DataOut(idx).EL_record       = EL_record;
    %DataOut(idx).ID_record       = ID_record;
    DataOut(idx).information      = information;
    %DataOut(idx).SV              = SV_record;
    DataOut(idx).IonoCoefficients = IonoCoefficients;
    DataOut(idx).Function         = 'SatelliteLiteV3';
    DataOut(idx).include          = include;

    DataOut(idx).RXa              = RXa;
    DataOut(idx).RXb              = RXb;
    DataOut(idx).SVg_e            = SV;
    %DataOut(idx).r_ea_e1         = mprofile.r_ea_e1;
    %DataOut(idx).r_eb_e1         = mprofile.r_eb_e1;
    
    
end

fprintf('... Finished ... \n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function delays L1C
% Satellite specific delays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delay=delays(Ir, Tr, Tha, Mpa, Rxa, Na, id, input, include, IonoCoeff)
    % L1C delay
    % Space delays - constant residual per epoch
    % Iono
    % Tropo
    % Thermal
    % Multipath
    % Receiver bias
    % Random
    % Integer ambiguity
    Ir_delay_a = Ir.delay(input,IonoCoeff)*include.Iono;
    Tr_delay_a = Tr.delay(input)          * include.Tropo; 
    Th_delay_a = Tha.delay(input)         * include.Th;  
    Mp_delay_a = Mpa(id).delay(input)     * include.MP;   
    Rx_delay_a = Rxa.delay                * include.Rx;   
    Rnd_delay_a= (randn * 0.1)./wgs84.c   * include.Random;
    
    
    Ir_delay_ca =  Ir.delay(input,IonoCoeff) * include.Iono;    
    Tr_delay_ca =  Tr.delay(input)           * include.Tropo;
    Th_delay_ca =  Tha.delay_ca(input)       * include.Th;
    Mp_delay_ca =  Mpa(id).delay_ca(input)   * include.MP;
    Rx_delay_ca =  Rxa.delay                 * include.Rx;
    Rnd_delay_ca= randn * 0.0                * include.Random;
    N_integer   =  Na(id).N                  * include.N; % cycles

    % Pseudorange delays
    delay.C1C.Ir  = Ir_delay_a;
    delay.C1C.Tr  = Tr_delay_a;
    delay.C1C.Th  = Th_delay_a;
    delay.C1C.Mp  = Mp_delay_a;
    delay.C1C.Rx  = Rx_delay_a;
    delay.C1C.Rnd = Rnd_delay_a;
    
    % Carrier phase delays
    delay.L1C.Ir  = Ir_delay_ca;
    delay.L1C.Tr  = Tr_delay_ca;
    delay.L1C.Th  = Th_delay_ca;
    delay.L1C.Mp  = Mp_delay_ca; 
    delay.L1C.Rx  = Rx_delay_ca;
    delay.L1C.Rnd = Rnd_delay_ca; 
    delay.L1C.N   = N_integer;
    
    % add include to delay
    delay.include = include;
end

function y=loadTrajectory(src, event)    
    % Iterated Weighted Least Squares : ECEF implementation
    % load Tracjectory that has a very specific format.
    % 1     [1x1 Signal]    time         October_Aircraft_ecef_ned/time         
    % 2     [1x1 Signal]    <Xe>         October_Aircraft_ecef_ned/xn xe xd     
    % 3     [1x1 Signal]    <U w>        October_Aircraft_ecef_ned/U w          
    % 4     [1x1 Signal]    ''           October_Aircraft_ecef_ned/phi theta psi
    % 5     [1x1 Signal]    <omega_ib_b> October_Aircraft_ecef_ned/omega_ib_b   
    % 6     [1x1 Signal]    <n>          October_Aircraft_ecef_ned/n            
    % 7     [1x1 Signal]    <windNED>    October_Aircraft_ecef_ned/windNED      
    % 8     [1x1 Signal]    CMD          October_Aircraft_ecef_ned/CMD          
    % 9     [1x1 Signal]    <f_ib>       October_Aircraft_ecef_ned/f_ib         
    % 10    [1x1 Signal]    <omega_ib>   October_Aircraft_ecef_ned/omega_ib     
    % 11    [1x1 Signal]    <lat lon>    October_Aircraft_ecef_ned/lat lon      
    % 12    [1x1 Signal]    <hd>         October_Aircraft_ecef_ned/hd           
    % 13    [1x1 Signal]    <Ve>         October_Aircraft_ecef_ned/Ve 
    %file = 'SIMDATA_AEROSONDE_100Hz_[f_ib_omega_ib]_TECS2_SND2_Long_09_08_19_wen.mat';
     
    load(config.file); 
    time  = simOut.time;
    dt    = diff(time);
    dt_avg= mean(dt);
    
    % TODO: specify receiver output rate and change subsampling factor
    % respectively
    odr   = 1;
    ssf   = round(odr/dt_avg); 
    
    time  = time(1:ssf:end);
    no_epochs = length(time);    
    
    mprofile.time  = time;                      % 1 second epoch
    mprofile.lat   = simOut.lat(1:ssf:end);     % Geodetic latitude  [radians]
    mprofile.lon   = simOut.lon(1:ssf:end);     % Geodetic longitude [radians]
    mprofile.hd    = simOut.hd(1:ssf:end);      % Ellipsoidal height [metres]
    mprofile.roll  = simOut.roll(1:ssf:end);    % Roll
    mprofile.ptch  = simOut.pitch(1:ssf:end);   % Pitch
    mprofile.yaw   = simOut.yaw(1:ssf:end);     % Yaw
    mprofile.v_ea_n= simOut.simOut.velNED(1:ssf:end).';% velocity in ned
    
    % geodetic2ecef
    a  = wgs84.A; % semimajor axis
    e2 = wgs84.f*(2-wgs84.f);
    
    % Get ECEF coordinates
    [x_ea_e, y_ea_e, z_ea_e]=ell2xyz(mprofile.lat,mprofile.lon,mprofile.hd,a,e2);  
    if size(x_ea_e, 2) == 1
        mprofile.r_ea_e1 = [x_ea_e, y_ea_e, z_ea_e];
    else
        mprofile.r_ea_e1 = [x_ea_e.', y_ea_e.', z_ea_e.'];
    end
    
    % ECEF velocity
    lat0 = simOut.lat(1);
    lon0 = simOut.lon(1);
    v_ea_e= zeros(size(mprofile.v_ea_n));
    for k=1:length(xN)
        v_ea_e(k,:) = (Rne([lat0, lon0])*mprofile.v_ea_n(k,:).').';
    end
    mprofile.v_ea_e  = v_ea_e;

    fprintf('Trajectory Loaded... \n');
    y=no_epochs;  
end

end







