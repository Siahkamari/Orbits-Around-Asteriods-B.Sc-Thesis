function dydt = gauss(~,C)

global mu k

a = C(1); e = C(2); RA = C(3);
incl = C(4); w = C(5); TA = C(6);

n = sqrt(mu/a^3);
r = a*(1-e^2)/(1+e*cos(TA));
u = w + TA;


% Propagating for a
da = -(4*a*e*k*cos(incl)*...
    (a*(-1 + e^2) + r + e*r*cos(TA))*sin(TA))/((-1 + e^2)*r);

% Propagating for e
de =-((2*(-1 + e^2)*k*cos(incl)*sin(TA))/(1 + e*cos(TA))) ;
% Propagating for RA
dRA = (2*k*r*sin(u)*(sin(u) + e*sin(w)))/(a*(-1 + e^2));

% Propagating for incl
di = (2*k*r*cos(u)*sin(incl)*(sin(u) + e*sin(w)))/(a*(-1 + e^2)) ;

% Propagating for w
dw = -(1/(a*e*(-1 + e^2)))*k*cos(incl)*...
    (-e*r*cos(2*u) + 2*a*(-1 + e^2)*cos(TA) +...
    e*r*cos(TA)^2 +... 
    e*(-2*a + 2*a*e^2 - r*sin(TA)^2 + 2*e*r*sin(u)*sin(w)));

% Propagating for TA
dTA = (a^2*sqrt(1 - e^2)*n)/r^2 + ...
     2*k*cos(incl)*...
     (((a*(-1 + e^2) - r)*sin(TA)^2)/(a*(-1 + e^2)) +...
     ((cos(TA)*(1 + e*cos(TA)))/e));

   % =================================================== %
   
% Propagating for a
da1 = (1/(sqrt(1 - e^2)*n))*...
  2*k^2*(a*(-1 + e^2)*cos(u)*sin(incl)^2*sin(u) +... 
   e*r*cos(u)^2*sin(TA) + e*r*cos(incl)^2*sin(u)^2*sin(TA));

% Propagating for e
de1 =(1/(a*n))*sqrt(1 - e^2)*k^2*r*...
 (-cos(u)*(cos(TA) + (e + cos(TA))/(1 + e*cos(TA)))*...
sin(incl)^2*sin(u) + (cos(u)^2 + cos(incl)^2*sin(u)^2)*sin(TA)) ;
% Propagating for RA
dRA1 = -((k^2*r^2*cos(incl)*sin(u)^2)/(a^2*sqrt(1 - e^2)*n))-0*k;

% Propagating for incl
di1 = -((k^2*r^2*cos(incl)*cos(u)*sin(incl)*sin(u))/(a^2*sqrt(1 - e^2)*n)) ;

% Propagating for w
dw1 = (1/(a^2*e*sqrt(1 - e^2)*n))*k^2*r*...
    (a*(-1 + e^2)*cos(u)^2*cos(TA) +... 
    cos(incl)^2*(e*r + a*(-1 + e^2)*cos(TA))*sin(u)^2 +...
    (a*(-1 + e^2) - r)*cos(u)*sin(incl)^2*sin(u)*sin(TA));

% Propagating for TA
dTA1 =0*(a^2*sqrt(1 - e^2)*n)/r^2 +...
   1/(a^2*e*sqrt(1 - e^2)*n)*k^2*r*...
   ((a - a*e^2 + r)*cos(u)*sin(incl)^2*sin(u)*sin(TA)-...
   a*(-1 + e^2)*cos(TA)*(cos(u)^2 + cos(incl)^2*sin(u)^2));

% =================================================== %


dydtco = [da;de;dRA;di;dw;dTA;dTA-(a^2*sqrt(1 - e^2)*n)/r^2];
dydtce = [da1;de1;dRA1;di1;dw1;dTA1;dTA1];
dydt = dydtco + dydtce;
end