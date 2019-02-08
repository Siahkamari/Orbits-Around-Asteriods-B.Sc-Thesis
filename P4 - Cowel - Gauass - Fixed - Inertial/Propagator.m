function y = Propagator(T,C)

global mu Jx Jy u

a = C(1);
e = C(2);
RA = C(3);
incl = C(4);
w = C(5);
TA = C(6);
n = sqrt(mu/a^3);
p = a*(1-e^2);
k = 0*3*n*(Jx+Jy)/(4*p^2);
c = 0*(Jx-Jy)/(Jx+Jy);


% Propagating for a
y(:,1) = a + 0*T;

% Propagating for e
y(:,2) = e + 0*T;

% Propagating for RA
y(:,3) = RA - k*cos(incl)*(1 - c*cos(2*RA))*T +...
    u^2*cos(incl)*(-2-3*e^2)/(4*n*(1-e^2)^0.5)*T-u*T;

% Propagating for incl
y(:,4) = incl +k*c*sin(incl)*sin(2*RA)*T;

% Propagating for w
y(:,5) = w + k/4*(3+5*cos(2*incl)-c*(-1+5*cos(2*incl))*cos(2*RA))*T;

% Propagating for M0
y(:,6) = (n+u*cos(incl))*T;


end