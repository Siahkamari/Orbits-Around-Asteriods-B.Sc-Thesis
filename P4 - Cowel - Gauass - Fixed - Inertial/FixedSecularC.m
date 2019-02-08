function dy = FixedSecularC(~,C)

global mu C20 C22 u

a = C(1); e = C(2); RA = C(3);
incl = C(4);

n = sqrt(mu/a^3);
p = a*(1-e^2);

k = 3*n/(8*p^2);

%% Propagating

da = 0;
de = 0;
dRA = 4*k*cos(incl)*(C20 + 2*C22*cos(2*RA)) - u;
di = 8*k*C22*sin(incl)*sin(2*RA);
dw = -k*(C20*(3+5*cos(2*incl))+ 2*C22*(-1+5*cos(2*incl))*cos(2*RA));
dTA = n + k*(C20*(1+3*cos(2*incl)) -12*C22*cos(2*RA)*sin(incl)^2) + u*cos(incl) ;

% =================================================== %

dy = [da;de;dRA;di;dw;dTA;u*cos(incl)];

end