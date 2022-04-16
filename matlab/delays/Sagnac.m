% Author    : Hery A Mwenegoha (C) UoN 2019
% version
%   current : 1.01    18/05/2019
%   first   : 1.00    01/02/2019
% Description
%   Sagnac modelling with receiver in either NED|Curvillinear or ECEF frame.
%   Satellite is always treated to be in *ECEF* frame.
% Functions
%   ned  - delay with receiver in NED  frame  [m]
%   ecef - delay with reciver  in ECEF frame  [m]
% syntax
%   y=ned(xned,    r_es_e)
%   y=ecef(r_ea_e, r_es_e)
% input
%   xned
%   Latitude  : mu receiever                  [rad]
%   Longitude : lambda receiver               [rad]
%   height    : ellipsoid height receiver     [m]
%   r_ea_e
%   Xe        : x-ecef receiver position      [m]
%   Ye        : y-ecef receiver position      [m]
%   Ze        : z-ecef receiver position      [m]
%   r_es_e
%   Xe        : x-ecef satellite position     [m]
%   Ye        : y-ecef satellite position     [m]
%   Ze        : z-ecef satellite position     [m]
% output
%   y         : time delay/advance due to sagnac effect [m]
classdef Sagnac < handle
    methods(Static)
        function y = metres(x, r_es_e)
            lat= x(1);
            lon= x(2);
            h  = x(3);
            Xs = r_es_e(1);
            Ys = r_es_e(2);
            Zs = r_es_e(3);
            y=-(wgs84.omega_ie*(Ys*cos(lat)*cos(lon)*(h + 6378137/(1 - (482380915665231*sin(lat)^2)/72057594037927936)^(1/2)) - Xs*cos(lat)*sin(lon)*(h + 6378137/(1 - (482380915665231*sin(lat)^2)/72057594037927936)^(1/2))))/wgs84.c_light;
        end
        
        function y = ned(x, r_es_e)
            lat= x(1);
            lon= x(2);
            h  = x(3);
            Xs = r_es_e(1);
            Ys = r_es_e(2);
            Zs = r_es_e(3);
            y=-(wgs84.omega_ie*(Ys*cos(lat)*cos(lon)*(h + 6378137/(1 - (482380915665231*sin(lat)^2)/72057594037927936)^(1/2)) - Xs*cos(lat)*sin(lon)*(h + 6378137/(1 - (482380915665231*sin(lat)^2)/72057594037927936)^(1/2))))/wgs84.c_light;
        end
        
        function y = ecef(r_ea_e, r_es_e)
            y=-wgs84.omega_ie/wgs84.c_light *...
                (r_es_e(2)*r_ea_e(1) - r_es_e(1)*r_ea_e(2));
        end
        
        function y = range(r_ea_e, r_es_e)
            % metres
            y=Sagnac.ecef(r_ea_e, r_es_e);
        end
        
        function y = rate(r_ea_e, v_ea_e, r_es_e, v_es_e)
            % metres/second
            y=-wgs84.omega_ie/wgs84.c_light *...
                ( v_es_e(2)*r_ea_e(1) +  r_es_e(2)*v_ea_e(1)...
                - v_es_e(1)*r_ea_e(2)- r_es_e(1)*v_ea_e(2));
        end
    end
end