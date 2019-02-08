%% Main *******(Havaset be density o vertexo shekl bashe ina bashe)
clear all
clc

global G density Face n
G = 6.67384e-11;
Mass = 5.97e24;
Re = 6378e3;
density = Mass/48e18;
duration = 86400;

 
%% Inputtig the shape
load('vertex_cube.mat');
load('face_cube.mat');
vertex = vertex*1e6;

n = numel(face(:,1));
Face = struct('Pi',[],'Pj',[],'Pk',[],...
    'eij',[],'ejk',[],'eki',[],'nf',[],'nij',[],'njk',[],'nki',[]);
Face(n) = struct('Pi',[],'Pj',[],'Pk',[],...
    'eij',[],'ejk',[],'eki',[],'nf',[],'nij',[],'njk',[],'nki',[]);


%% Getting the degree from user
for f = 1:n
    
    j1 = face(f,2);
    j2 = face(f,3);
    j3 = face(f,4);
    Pi = [vertex(j1,2),vertex(j1,3),vertex(j1,4)];
    Pj = [vertex(j2,2),vertex(j2,3),vertex(j2,4)];
    Pk = [vertex(j3,2),vertex(j3,3),vertex(j3,4)];
    
    eij = norm(Pj-Pi);
    ejk = norm(Pk-Pj);
    eki = norm(Pi-Pk);
    
    nf = cross((Pj-Pi),(Pk-Pi));
    nf = nf/norm(nf);
    
    nij = cross((Pj-Pi),nf);
    nij = nij/norm(nij);
    
    njk = cross((Pk-Pj),nf);
    njk = njk/norm(njk);
    
    nki = cross((Pi-Pk),nf);
    nki = nki/norm(nki);
    
    Face(f) = struct('Pi',Pi,'Pj',Pj,'Pk',Pk,...
        'eij',eij,'ejk',ejk,'eki',eki,'nf',nf,'nij',nij,'njk',njk,'nki',nki);
end
clear face vertex j1 j2 j3 Pi Pj Pk eij ejk eki nf nij njk nki f


%% Getting the Trinomials for shperical coefficients

X0 = [7000,0,0,0,2,8]*1e3;
options=odeset('Reltol',10^-5);
[T,X]=ode45(@X_dot,[0 duration],X0,options);
R=[X(:,1),X(:,2),X(:,3)];
V=[X(:,4),X(:,5),X(:,6)];


%% Printimng the body
figure(1)
plot3(R(:,1),R(:,2),R(:,3))
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
grid on;
hold all

a = 1000000;
b = 2000000;
c = 3000000;
plotcube(a,b,c)

axis equal
grid on
