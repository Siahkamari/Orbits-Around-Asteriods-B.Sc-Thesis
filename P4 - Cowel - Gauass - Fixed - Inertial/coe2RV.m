function [r, v] = coe2RV(coe)

global mu u

a = coe(1); e = coe(2); RA = coe(3);
incl = coe(4); w = coe(5); TA = coe(6);
rp = (a*(1-e^2)) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (sqrt(mu/(a*(1-e^2)))) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

R3_W = [cos(RA) sin(RA) 0
        -sin(RA) cos(RA) 0
        0 0 1];

R1_incl = [1 0 0
        0 cos(incl) sin(incl)
        0 -sin(incl) cos(incl)];

R3_w = [cos(w) sin(w) 0
        -sin(w) cos(w) 0
        0 0 1];
    
% transmition matrix from perifocal into inretial
Q_PQW2GCRF = (R3_w*R1_incl*R3_W)';

r = Q_PQW2GCRF*rp;
v = Q_PQW2GCRF*vp;

r = r';
v = v';
end