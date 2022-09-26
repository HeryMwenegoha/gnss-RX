% Hery A Mwenegoha copyright (c) 2020

function y = quat2euler(q)
% Quaternion to Euler
q0    = q(1);
q1    = q(2);
q2    = q(3);
q3    = q(4);
phi   = atan2(2*q0*q1 + 2*q2*q3, q0^2 - q1^2 - q2^2 + q3^2);
theta = asin(2*q0*q2 - 2*q1*q3);
psi   = atan2(2*q0*q3 + 2*q1*q2, q0^2 + q1^2 - q2^2 - q3^2);
y     = [phi;theta;psi];
end