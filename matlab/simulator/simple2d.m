% Hery A Mwenegoha copyright (c) 2020

% Initial position
lat0 = deg2rad(52.9519817);
lon0 = deg2rad(-1.1884625);
h0   = 0;

% ECEF at navOrigin
[X0, Y0, Z0] = geodetic2ecef(lat0, lon0, h0, 'wgs84');

% Trajectory in local nav
north = [0, 100, 100,    0,    0, 100];
east  = [0,   0, 100,  100,  200, 200];
down  = [0,   0,   0,    0,    0,   0];
t = linspace(0,60, length(north));

% time
time   = t(1):0.01:t(end);

% Interpolation
posNorth = spline(t, north, time);
posEast  = spline(t, east, time);
posDown  = spline(t, down, time);

% Positions in ECEF
lat = [];
lon = [];
hd  = [];
for i=1:length(posNorth)
Rne    = Ren([lat0, lon0]).';
dEcef  = Rne * [posNorth(i); posEast(i); posDown(i)];
ecefPos= [X0; Y0; Z0] + dEcef;
[lat(end+1),lon(end+1), hd(end+1)] = ecef2geodetic(ecefPos(1), ecefPos(2), ecefPos(3), 'wgs84');
end

% velocities
vN = diff(posNorth)./diff(time); 
vE = diff(posEast)./diff(time); 
vD = diff(posDown)./diff(time); 

vN = [vN(1) vN];
vE = [vE(1) vE];
vD = [vD(1) vD];

roll = zeros(1, length(vE));
pitch= zeros(1, length(vE));
yaw  = atan2(vE, vN);

str.time  = time;         % 1 second epoch
str.lat   = lat;          % Geodetic latitude   [radians]
str.lon   = lon;          % Geodetic longitude  [radians]
str.hd    = hd;           % Ellipsoidal height  [metres]
str.roll  = roll;         % Roll                [radians]
str.pitch = pitch;        % Pitch               [radians]
str.yaw   = yaw;          % Yaw                 [radians]
str.velNED= [vN; vE; vD]; % velocity in ned     [metres-per-second]

%save('simulator/simple2d-03-04.mat','str');

figure('Position', [440 132 569 665]);
subplot(3,1,1);
plot(posEast,posNorth, 'LineWidth', 2); hold on;
plot(east,   north,'*', 'MarkerSize', 15); hold off;
xlabel('East [m]');
ylabel('North [m]');
grid minor;

subplot(3,1,2);
plot(time, vN, 'LineWidth', 2); hold on;
plot(time, vE, 'LineWidth', 2); hold off;
xlabel('time [s]');
grid minor;
ylabel('Velocity [m/s]');
legend('vN', 'vE');

subplot(3,1,3);
plot(time, rad2deg(roll),  'LineWidth', 2); hold on;
plot(time, rad2deg(pitch), 'LineWidth', 2); 
plot(time, rad2deg(yaw),   'LineWidth', 2); hold off;
xlabel('time [s]');
grid minor;
ylabel('\phi, \theta, \psi [deg]');

