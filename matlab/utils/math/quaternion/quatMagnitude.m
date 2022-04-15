% Hery A Mwenegoha copyright (c) 2020

function y=quatMagnitude(q)
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);
y = sqrt(q0.^2 + q1.^2 + q2.^2 + q3.^2);
end