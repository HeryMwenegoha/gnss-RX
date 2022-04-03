% Hery A Mwenegoha copyright (c) 2020

function [mprofile, y]=loadTrajectory(trajectoryFilename)
% loads a trajectory with str-structure:
% str: time
%    : lat
%    : lon
%    : hd
%    : roll
%    : pitch
%    : yaw
%    : v_ea_n : NED  position of the receiver
%    : r_ea_e : ECEF position of the receiver
%    : v_ea_e : ECEF velocity of the receiver
load(trajectoryFilename);
time  = str.time;
dt    = diff(time);
dt_avg= mean(dt);

% TODO: specify receiver output rate and change subsampling factor
% respectively
odr   = 1;
ssf   = round(odr/dt_avg);

time  = time(1:ssf:end);
no_epochs = length(time);

mprofile.time  = time;                      % 1 second epoch
mprofile.lat   = str.lat(1:ssf:end);     % Geodetic latitude  [radians]
mprofile.lon   = str.lon(1:ssf:end);     % Geodetic longitude [radians]
mprofile.hd    = str.hd(1:ssf:end);      % Ellipsoidal height [metres]
mprofile.roll  = str.roll(1:ssf:end);    % Roll
mprofile.ptch  = str.pitch(1:ssf:end);   % Pitch
mprofile.yaw   = str.yaw(1:ssf:end);     % Yaw
mprofile.v_ea_n= str.velNED(:,1:ssf:end);% velocity in ned

% Get ECEF coordinates
% x_ea_e : receiver a x-ecef-pos
% y_ea_e : receiver a y-ecef-pos
% z_ea_e : receiver a z-ecef-pos
[x_ea_e, y_ea_e, z_ea_e]=geodetic2ecef(mprofile.lat,mprofile.lon,mprofile.hd, 'wgs84');
mprofile.r_ea_e = [x_ea_e; y_ea_e; z_ea_e];


% ECEF velocity
lat0 = str.lat(1);
lon0 = str.lon(1);
v_ea_e= zeros(size(mprofile.v_ea_n));
for k=1:length(mprofile.v_ea_n)
    v_ea_e(:,k) = Rne([lat0, lon0])*mprofile.v_ea_n(:,k);
end
mprofile.v_ea_e  = v_ea_e;

fprintf('Trajectory Loaded... \n');
y=no_epochs;
end