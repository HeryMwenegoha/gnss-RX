% Hery A Mwenegoha copyright (c) 2020

function [gnssObj,rawxSoln] = gnss_rx_update(gnssObj, current_SOW, posEcef, velEcef, Rpy)
% This function is called every navEpoch to compute the raw GNSS
% observables. 
% Inputs:
%   gnssObj - the gnssObj structure containing information for the
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

% Current epoch g_time
t=current_SOW;

% Get the current receiver position from the Navigation engine
r_ea_e1   = posEcef; % 3x1
v_ea_e1   = velEcef; % 3x1
att.roll  = Rpy(1);  % 1x1
att.ptch  = Rpy(2);  % 1x1
att.yaw   = Rpy(3);  % 1x1

% TODO: include GALILEO, GLONASS, etc
% Get the space vehicles
svGps = gnssObj.svGps;

% TODO: we definitely want a way of checking this per navEpoch. But the
% Nav Engine epoch (incoming pose solution) and gnss-RX epochs are
% different for our case.
% Number of Satellites in SV struct
nSats  = length(svGps);

% TODO: Get better DOY from code below
% HACK-DOY for the epoch
DOY = gnssObj.config.DOY;
% Get the DOY
% sow_lapsed= t-time_flght_start.g_sow;
% g_t2 = time_flght_start + seconds(sow_lapsed);
% DOY = g_t2.g_doy;

% calculate satellite positions at time of transmission - t
for nPRN=1:nSats
    if isempty(svGps(nPRN).ID) == true
        continue;
    end
    
    id = svGps(nPRN).ID;
    
    % Get classes
    class.config    = gnssObj.config;
    class.include   = gnssObj.include;
    class.IonoCoeff = gnssObj.IonoCoeff;
    class.Ir        = gnssObj.Ir;
    class.Tr        = gnssObj.Tr;
    class.Mp        = gnssObj.Mpa(id);
    class.Th        = gnssObj.Tha;
    class.Rx        = gnssObj.Rxa;
    class.Nn        = gnssObj.Na(id);
    class.r_e       = r_ea_e1;
    class.v_e       = v_ea_e1;
    class.att       = att;
    
    % Compute measurements
    [gT_r,C1C,L1C,D1C,LLI]=sat.PR(svGps, id,t,class);
    
    % Store measurements if OK
    if ~isnan(C1C)
        % Receiver data holder
        % -- LLI - flag
        rawxSoln.gTr = gT_r;
        rawxSoln.SOW = gT_r*1;       % can be presented ok
        rawxSoln.DOY = DOY;
        rawxSoln.svg(id).PRN = svGps(id).ID;
        rawxSoln.svg(id).C1C = C1C;
        rawxSoln.svg(id).L1C = L1C;
        rawxSoln.svg(id).D1C = D1C;
        rawxSoln.svg(id).S1C = [];
        rawxSoln.svg(id).LLI = LLI;
        rawxSoln.bias_m      = gnssObj.Rxa.delay*wgs84.c;
        rawxSoln.drift_mps   = gnssObj.Rxa.ddelay*wgs84.c;
        rawxSoln.pos_ecef    = r_ea_e1;
        rawxSoln.vel_ecef    = v_ea_e1;        
    end
    
end

% RX-A CLOCK MODEL
gnssObj.Rxa.common;

% IONO-GM Process
gnssObj.Ir.common;

% TROPO-GM Process
gnssObj.Tr.common;
end