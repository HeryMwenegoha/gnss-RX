%% Call constructors

% Other constructors

% Initialise receiver gnssObj
gnssObj = gnss_rx_constructor();

% Load a fake trajectory
[mprofile, maxSeconds] = loadTrajectory(gnssObj.config.userFile);

% Get current SOW
tStart = gnssObj.config.start_g_time_s.g_sow;
SOW    = tStart:tStart+60;

%% Main loop
fprintf('Start soln \n');
tic;
for navEpoch = 1:length(SOW)
    % current time sec
    current_time_s = SOW(navEpoch);
    
    % Get the current navSolution that will be used to solve for the raw
    % GNSS solution.
    if exist('mprofile', 'var')
        lat = mprofile.lat(navEpoch);
        lon = mprofile.lon(navEpoch);
        hd  = mprofile.hd(navEpoch);
        velNED = mprofile.v_ea_n(:,navEpoch);
        Rpy = [mprofile.roll(navEpoch);
            mprofile.ptch(navEpoch);
            mprofile.yaw(navEpoch)];
    else
        lat = deg2rad(52.9519816);
        lon = deg2rad(-1.1907585);
        hd  = 100;
        velNED= [0;0;0];
        Rpy = [0;0;0];
    end
    
    % Get ECEF coordinates from simulated llh
    [xEcef, yEcef, zEcef] = geodetic2ecef(lat, lon, hd, 'wgs84');
    
    % Pack 
    posEcef = [xEcef;yEcef;zEcef];
    velEcef = Ren(lat, lon).'*velNED;
    
    % Pass variables into update function
    [gnssObj,rawxSoln] = gnss_rx_update(gnssObj, current_time_s, posEcef, velEcef, Rpy);
end
toc;
fprintf('End soln \n');