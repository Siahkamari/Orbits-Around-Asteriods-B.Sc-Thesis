function dy = FixedSecular(~,C)

global mu Jx Jy u


a = C(1); e = C(2); RA = C(3);
incl = C(4);

n = sqrt(mu/a^3);
p = a*(1-e^2);

k = 3*n*(Jx+Jy)/(4*p^2);
c = (Jx-Jy)/(Jx+Jy);


% Propagating for a
da = 0;

% Propagating for e
de = 0;
% Propagating for RA
dRA = - k*cos(incl)*(1 - c*cos(2*RA)) - u;

% Propagating for incl
di = k*c*sin(incl)*sin(2*RA);

% Propagating for w
dw = k/4*(3+5*cos(2*incl)-c*(-1+5*cos(2*incl))*cos(2*RA));

% Propagating for TA
dTA = n ; %***** bayad be ezafeye chizi shavad*******;

% =================================================== %

dy = [da;de;dRA;di;dw;dTA];

end