%% Main
clear
clc

str = {'Cowel','InertialSecular','FixedSecular','Gauss(no perturbation)'};
[method,ok] = listdlg('PromptString','Select a Method:',...
    'SelectionMode','single',...
    'ListSize',[250,120],...
    'ListString',str);

if ok == 0
    clear
    break
end

%% Constants
global mu Jx Jy u k duration C20 C22
% mu = 5.60616e6;
% u = 0.000001;
% Jx = 2e8;
% Jy = 0.35e8;
% C22 = 0.25*(Jx - Jy);
% C20 = -0.5*(Jy + Jx);

%% Fixed orbit
Jx = 600e3;
Jy = -200e3;
C22 = 0.25*(Jx - Jy);
C20 = -0.5*(Jy + Jx);
mu = 5e12*6.674e-11;
u = 5e-6;
% X0 = [0,150000,0,2.088-0*2e5*u,0,4.8];

% mu = 398600.4415;
% u = 366/365*(2*pi)/(86400);
% 433 Eros
% mu = 5e-4;
% C20 = -26.755;
% C22 = 12.725;
% u = 3.3118e-4;
% Itokawa
% mu = 4.01e-9;
% C20 = -6.6133e-3;
% C22 = 4.7972e-3;
% u = 1.4386e-4;
% Jx = -C20 + 2*C22;
% Jy = -C20 - 2*C22;

k = u;
duration = 1000000/4;

% X0 = [100000,120000,55,1.5,2,4.8];

coe_0 = [3434,0.04,pi/2,132.5*pi/180,0,0];
% coe_0(3:6) = pi/180*coe_0(3:6);
[R,V] = coeFixed2RV(coe_0,0);
X0 = [R,V];
coei_0 = RV2coe(R,V);
options=odeset('Reltol',10^-10,'abstol',10^-10);

%% Choosing Method
if method == 1
    %% Cowel
    [T,X]=ode45(@X_dot,0:500:duration,X0,options);
    R=[X(:,1),X(:,2),X(:,3)];
    V=[X(:,4),X(:,5),X(:,6)];
    color = 'b';
    
    n = numel(T);
    coe = zeros(n,6);
    
    for i = 1:n
        coe(i,:) = RV2coe(R(i,:),V(i,:));
    end
    coei = coe;
  
    options=odeset('Reltol',10^-13,'abstol',10^-16);
elseif method==2
    
    %% Inertial Secular
%     coe_0 = RV2coe(X0(1:3),X0(4:6));
coei_0(1) = 3480;
coei_0(2) = 0.03;
coei_0(4) = 131*pi/180;
    [T,coei]=ode45(@InertialSecular,[0 duration],coei_0,options);
    color = 'g';
    n = numel(T);
    R = zeros(n,3);
    V = zeros(n,3);
    for i = 1:n
        [R(i,:),V(i,:)] = coe2RV(coei(i,:));
    end
    
elseif method == 4
    %% Fixed Frame (guass)
%     X0 = X0+[0,0,0,X0(2)*u,-X0(1)*u,0];
%     coe_0 = RV2coe(X0(1:3),X0(4:6));
    coe_0(7) = coe_0(6);
    [T,coe]=ode45(@gauss,0:500:duration,coe_0,options);
    n = numel(T);
    color = 'b';
    coei = zeros(n,6);
    R = zeros(n,3);
    V = zeros(n,3);
    
    for i = 1:n
        [R(i,:),V(i,:)] = coeFixed2RV(coe(i,:),T(i));
    end
    
    for i = 1:n
        coei(i,:) = RV2coe(R(i,:),V(i,:));
    end
    
else
    %% Fixed Frame (guass and Secular can be used)
%     X0 = X0+[0,0,0,X0(2)*u,-X0(1)*u,0];
%     coe_0 = RV2coe(X0(1:3),X0(4:6));
    coe_0(7) = coe_0(6);
    coe_0(1) = 3500;
    [T,coe]=ode45(@FixedSecularC,0:500:duration,coe_0,options);
    n = numel(T);
    color = 'r';
    coei = zeros(n,6);
    R = zeros(n,3);
    V = zeros(n,3);
    
    for i = 1:n
        [R(i,:),V(i,:)] = coeFixed2RV(coe(i,:),T(i));
    end
    
    for i = 1:n
        coei(i,:) = RV2coe(R(i,:),V(i,:));
    end
end
% %% BE Carefullllllllllllllllllllllllllll
% Rf=R;
% for i=1:n
%     t=T(i);
%     Q_Fixed2Inertial = [cos(u*t) -sin(u*t)  0;sin(u*t) cos(u*t)  0;0    0   1];
%     Rf(i,:)= Q_Fixed2Inertial'*R(i,:)';
% end
% R=Rf;
%%
% Ploting orbit
figure(1)
plot3(R(:,1),R(:,2),R(:,3),color);
xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
grid on;
hold all
%[x,y,z] = sphere;
%mesh(6378*x,6378*,6358*z);
%axis tight;

% Plotting Inertial Elements
figure(2)
Ploting_coe2(coei,T,color)
% end
% legend('Cowel method','Slowly rotating asteroid method','Fixed-body method')