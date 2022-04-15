% Hery A Mwenegoha copyright (c) 2020

function y = quat2Rbn(q)
% quat2Rbn : qbn2Rbn
% Quaternion to DCM-NED
% No difference to quat2Rbe just notation
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);
y  = [q0.^2+q1.^2-q2.^2-q3.^2  2*(q1*q2 - q3*q0)        2*(q1*q3 + q2*q0);
    2*(q1*q2 + q3*q0)         q0.^2-q1.^2+q2.^2-q3.^2  2*(q2*q3 - q1*q0);
    2*(q1*q3 - q2*q0)         2*(q2*q3 + q1*q0)        q0.^2-q1.^2-q2.^2+q3.^2];
end