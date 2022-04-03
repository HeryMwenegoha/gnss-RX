% Hery A Mwenegoha copyright (c) 2020

function  y = wrap_360(deg)
res     = mod(deg, 360);
ip      = res < 0;
res(ip) = res(ip)+ 360;
y       = res;
end