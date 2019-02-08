function dXdt = X_dot(t,X)
t
global mu Jx Jy u

% Some useful variables
R = X(1:3);
V = X(4:6);
x = cos(u*t)*X(1) + sin(u*t)*X(2);
y = -sin(u*t)*X(1) + cos(u*t)*X(2);
% x = X(1);
% y = X(2);
z = X(3);
r = norm(R);

% Defining the perturbed potential
Ux = 3*mu*x*(3*Jx+Jy)/(2*r^5) -15*mu*x*(Jx*x^2+Jy*y^2)/(2*r^7);
Uy = 3*mu*y*(Jx+3*Jy)/(2*r^5) -15*mu*y*(Jx*x^2+Jy*y^2)/(2*r^7);
Uz = 3*mu*z*(Jx+Jy)/(2*r^5) -15*mu*z*(Jx*x^2+Jy*y^2)/(2*r^7);

% The equation of motion in rotating frame

% apx = 2*u*V(2) + u^2*x + Ux;
% apy = -2*u*V(1) + u^2*y + Uy;
% apz = Uz;

%  
apx = cos(u*t)*Ux - sin(u*t)*Uy;
apy = sin(u*t)*Ux + cos(u*t)*Uy;
apz = Uz;


% Vectorial acceleration
a2body = -mu*(R/r^3);
ap = [apx;apy;apz];

% Passing to output
dXdt = [V;a2body+ap];

end