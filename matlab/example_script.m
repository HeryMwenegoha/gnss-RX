%% Call constructors

% Output data rate. Included in every delay-model constructor
% inside Gnss.m (RXbias/Iono/Tropo/RXthermal/Multipath) - set to > 1 Hz
% here specifically to exercise the high-rate path.
odr_Hz = 10;

% Initialise receiver
Rcv = Gnss('odr_Hz', odr_Hz);

% Load a simulated trajectory, resampled to the receiver's own odr_Hz.
motionFile = fullfile('simulator', 'simple2d-03-04.mat');
[mprofile, maxSeconds] = loadTrajectory(motionFile, odr_Hz);

% get current SOW
tStart = Rcv.gpsStartTime_s.g_sow;
SOW    = tStart + mprofile.time;

% numsats and numEpochs
nSats   = length(Rcv.svGps);
nEph    = length(SOW);

% simple logging items
prMat   = nan(nEph, nSats);
doppMat = nan(nEph, nSats);
carrMat = nan(nEph, nSats);
cn0     = nan(nEph, nSats);
prnList = [];

%% Main loop
fprintf('Start soln \n');
tic;
for iEph = 1:length(SOW)
    % current time sec
    current_time_s = SOW(iEph);

    % use provided motionFile otherwise use fixed location
    if exist('mprofile', 'var')
        lat    = mprofile.lat(iEph);
        lon    = mprofile.lon(iEph);
        hd     = mprofile.hd(iEph);
        velNED = mprofile.v_ea_n(:,iEph);
        Rpy    =[mprofile.roll(iEph);
                 mprofile.ptch(iEph);
                 mprofile.yaw(iEph)];
    else
        lat    = deg2rad(52.9519816);
        lon    = deg2rad(-1.1907585);
        hd     = 100;
        velNED = [0;0;0];
        Rpy    = [0;0;0];
    end

    % get ECEF coordinates from llh
    [xEcef,...
     yEcef,...
     zEcef] = geodetic2ecef(lat, lon, hd, 'wgs84');
    posEcef = [xEcef;yEcef;zEcef];
    velEcef = Ren(lat, lon).'*velNED;

    % pass variables into the update function
    rawxSoln = Rcv.update(current_time_s, posEcef, velEcef, Rpy);

    % get observables
    C1C = {rawxSoln.svg.C1C};
    D1C = {rawxSoln.svg.D1C};
    L1C = {rawxSoln.svg.L1C};
    Cn0 = {rawxSoln.svg.cn0};
    PRN = {rawxSoln.svg.PRN};

    % get masking field (choose non-empty)
    ip  = cellfun(@(c)~isempty(c), C1C);

    if ~any(ip)
        fprintf(2, 'no valid observables this epoch : %.2f \n', current_time_s-SOW(1));
        continue;
    end

    % log data
    prns              = [PRN{ip}];
    prMat(iEph,   ip) = [C1C{ip}];
    doppMat(iEph, ip) = [D1C{ip}];
    carrMat(iEph, ip) = [L1C{ip}];
    cn0(iEph,     ip) = [Cn0{ip}];

    if isempty(prnList)
        prnList = prns;
    else
        ipNotInList = ~ismember(prns, prnList);
        if any(ipNotInList)
            prnList = [prnList, prns(ipNotInList)];
        end
    end
end
toc;
fprintf('End soln \n');
fprintf('odr_Hz requested: %u , epoch spacing achieved: %.3f s\n', ...
    odr_Hz, median(diff(SOW)));


%%
figure('Name', 'Pr-Dopp-Carr-Consistency');

allPRNs  = [12, 14];

dt       =  median(diff(SOW));
diffPr   =  diff(prMat(:,allPRNs))./dt.*(1/0.1903);
diffPr   = -[diffPr(1,:); diffPr];
diffCarr =   diff(carrMat(:,allPRNs))./dt;
diffCarr = -[diffCarr(1,:);diffCarr];
dopp     =  doppMat(:,  allPRNs);

tiledlayout(2,2);
epoch   = SOW - SOW(1);
for isp = 1:length(allPRNs)
    nexttile
    hold on;
    plot(epoch, diffPr(:,   isp),   'LineWidth', 2, 'DisplayName', 'Pr');
    plot(epoch, dopp(:,  isp).', 'LineWidth', 2, 'DisplayName', 'f_D');
    plot(epoch, diffCarr(:, isp),   'LineWidth', 2, 'DisplayName', '\Phi')
    grid minor;
    xlabel('time [s]');
    legend;
    hold off;
    ylabel('cycles/s');
    xlabel('SOW [s]');
end
nexttile(3, [1,2])
plot(epoch, cn0(:, allPRNs), 'LineWidth',2)
xlabel('time [s]');
ylabel('cn0 [dB-Hz]');
grid minor;
legend(compose('G%02u', allPRNs))