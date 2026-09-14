% % epoch format for 2 receivers mounted on the aircraft
% fmt = struct('PRN',cell(1,32),'C1C',[],'L1C',[],'D1C',[],'S1C',[],'LLI',[]);
% RXa = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);
% RXb = struct('SOW',cell(1,600),'DOY',[],'svg',fmt,'svr',fmt);

% Initialise receiver gnssObj
% Both constructors draw randomly (iono-coefficient perturbation, and
% whatever internal state RXbias/Iono/Tropo/Multipath/Nr seed themselves
% with) - reseed between them so the two receivers start from identical
% initial conditions. Without this, the per-epoch sync in the loop below
% only matches what happens *inside* update(); it can't undo a mismatch
% baked in at construction, and since these are integrated/GM processes
% that mismatch persists for the whole run instead of averaging out.
odr_Hz  = 5;
rng(1);
gnssObj = gnss_rx_constructor('odr_Hz',odr_Hz);
rng(1);
G       = Gnss('odr_Hz',odr_Hz);

% Load a trajectory
[mprofile, maxSeconds] = loadTrajectory(gnssObj.config.userFile, odr_Hz);

% Get current SOW
tStart = gnssObj.config.start_g_time_s.g_sow;
SOW    = tStart+mprofile.time;

% pack to array
nSats   = length(gnssObj.svGps);
nEph    = length(SOW);

% prMat   = nan(nSats, length(SOW));
% doppMat = nan(nSats, length(SOW));
% carrMat = nan(nSats, length(SOW));
% cn0     = nan(nSats, length(SOW));

prMatOld   = nan(nEph, nSats);
doppMatOld = nan(nEph, nSats);
carrMatOld = nan(nEph, nSats);
cn0Old     = nan(nEph, nSats);

prMatNew   = nan(nEph, nSats);
doppMatNew = nan(nEph, nSats);
carrMatNew = nan(nEph, nSats);
cn0New     = nan(nEph, nSats);

prnList = [];

fprintf('Start soln \n');
tic;
for iEph = 1:nEph
    % current time sec
    current_time_s = SOW(iEph);
    
    % Get the current navSolution that will be used to solve for the raw
    % GNSS solution.
    if exist('mprofile', 'var')
        lat    = mprofile.lat(iEph);
        lon    = mprofile.lon(iEph);
        hd     = mprofile.hd(iEph);
        velNED = mprofile.v_ea_n(:,iEph);
        Rpy    = [mprofile.roll(iEph);
                  mprofile.ptch(iEph);
                  mprofile.yaw(iEph)];
    else
        lat    = deg2rad(52.9519816);
        lon    = deg2rad(-1.1907585);
        hd     = 100;
        velNED = [0;0;0];
        Rpy    = [0;0;0];
    end
    
    % Get ECEF coordinates from simulated llh
    [xEcef, yEcef, zEcef] = geodetic2ecef(lat, lon, hd, 'wgs84');
    %[xEcef, yEcef, zEcef] = ell2xyz(lat, lon, hd);
    
    % Pack 
    posEcef = [xEcef;yEcef;zEcef];
    velEcef = Ren(lat, lon).'*velNED;
    
    % Pass variables into both update functions, driven off the SAME
    % random draws: snapshot the generator before the old (struct-based)
    % scheme consumes it, then restore that snapshot before the new
    % (Gnss.m) scheme runs. Without this, RxaClk/Ir/Tr in the two
    % receivers are independent Gauss-Markov realizations that drift
    % apart over the run - the comparison below would be meaningless.
    rngState = rng;
    [gnssObj, rawxSolnOld] = gnss_rx_update(gnssObj, current_time_s, posEcef, velEcef, Rpy);

    rng(rngState);
    rawxSolnNew = G.update(current_time_s, posEcef, velEcef, Rpy);

    for iSat=1:nSats
        id = G.svGps(iSat).ID;

        if isempty(id) == true || isempty(rawxSolnOld.svg(id).C1C) || isempty(rawxSolnNew.svg(id).C1C)
            continue;
        end

        if ~any(prnList == id)
            prnList(end+1) = id;
        end

        prMatOld(iEph,   iSat) = rawxSolnOld.svg(id).C1C;
        doppMatOld(iEph, iSat) = rawxSolnOld.svg(id).D1C;
        carrMatOld(iEph, iSat) = rawxSolnOld.svg(id).L1C;
        cn0Old(iEph,     iSat) = rawxSolnOld.svg(id).cn0;

        prMatNew(iEph,   iSat) = rawxSolnNew.svg(id).C1C;
        doppMatNew(iEph, iSat) = rawxSolnNew.svg(id).D1C;
        carrMatNew(iEph, iSat) = rawxSolnNew.svg(id).L1C;
        cn0New(iEph,     iSat) = rawxSolnNew.svg(id).cn0;
    end
end
toc;
fprintf('End soln \n');

% Sanity: with the RNG synced above, Gnss.m should reproduce
% gnss_rx_constructor + gnss_rx_update to floating-point precision, not
% just "same order as the noise".
fprintf('max |old-new| C1C : %.3e m\n',   max(abs(prMatOld(:)   - prMatNew(:)),   [], 'omitnan'));
fprintf('max |old-new| L1C : %.3e cyc\n', max(abs(carrMatOld(:) - carrMatNew(:)), [], 'omitnan'));
fprintf('max |old-new| D1C : %.3e Hz\n',  max(abs(doppMatOld(:) - doppMatNew(:)), [], 'omitnan'));

%%
figure('Name', 'Pr-Dopp-Carr-Consistency');

allPRNs  = [12, 14];

dt       =  median(diff(SOW));
diffPr   =  diff(prMatNew(:,allPRNs))./dt.*(1/0.1903);
diffPr   = -[diffPr(1,:); diffPr];
diffCarr =   diff(carrMatNew(:,allPRNs))./dt;
diffCarr = -[diffCarr(1,:);diffCarr];
dopp     =  doppMatNew(:,  allPRNs);

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
plot(epoch, cn0New(:, allPRNs), 'LineWidth',2)
xlabel('time [s]');
ylabel('cn0 [dB-Hz]');
grid minor;
legend(compose('G%02u', allPRNs))
%% Old (gnss_rx_constructor/update) vs New (Gnss.m) consistency check
% With the RNG synced per-epoch above, these three panels should sit at
% ~0 (floating-point level), not at the metre/Hz-scale divergence you get
% when the two schemes are run without syncing the generator.
figure('Name', 'Old-vs-New-Consistency');
tiledlayout(1,3);

nexttile;
plot(epoch, prMatOld(:,allPRNs) - prMatNew(:,allPRNs), 'LineWidth', 2);
grid minor; xlabel('SOW [s]'); ylabel('\Delta C1C [m]');
title('Pseudorange: old - new'); legend(string(allPRNs));

nexttile;
plot(epoch, doppMatOld(:,allPRNs) - doppMatNew(:,allPRNs), 'LineWidth', 2);
grid minor; xlabel('SOW [s]'); ylabel('\Delta D1C [Hz]');
title('Doppler: old - new'); legend(string(allPRNs));

nexttile;
plot(epoch, carrMatOld(:,allPRNs) - carrMatNew(:,allPRNs), 'LineWidth', 2);
grid minor; xlabel('SOW [s]'); ylabel('\Delta L1C [cyc]');
title('Carrier-phase: old - new'); legend(string(allPRNs));