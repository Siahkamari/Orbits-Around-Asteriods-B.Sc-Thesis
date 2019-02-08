%================================ Main ===================================%

% This Programme calculates the Normalized spherical harmonics coefficients
% for the potential of a constant density polyhedron and Plots the shape

% Inputs : 
% vertex = the vertices of the polyhedra in matrix form
% face = the number of the vertices of each triangular plat in matrix form
% n = biggest degree for the coefficients to be calculated (input by user)
% R = Normalizing radius

% Outputs:
% The Volume of the body
% Spherical harmonics coefficients

% If we want Unnormalized coefficients we can Just change the k
% coefficients in the function Coefficients

clear
clc

%========================= Inputtig the shape ============================%
load('vertex_eros.mat');
load('face_eros.mat');
nv = numel(vertex(:,1));
nf = numel(face(:,1));


%========================= Normalizing radius ============================%
R = 16;


%==================== Getting the degree from user =======================%
fprintf('Input the degree (The nomalized coefficients will')
fprintf(' be computed from 0 to this degree and order) : \n')
n = input('');
clc


volume = 0;


for f = 1:nf
    
    j1 = face(f,2);
    j2 = face(f,3);
    j3 = face(f,4);
    x = [vertex(j1,2),vertex(j2,2),vertex(j3,2)];
    y = [vertex(j1,3),vertex(j2,3),vertex(j3,3)];
    z = [vertex(j1,4),vertex(j2,4),vertex(j3,4)];
    
    J = [x;y;z];
    
    volume = volume + det(J)/6;
end


%========== Getting the Trinomials for shperical coefficients ============%
[C,S] = Coefficients(n);
c = zeros(n+1,n+1);
s = zeros(n+1,n+1);


%=== Transfering coefficient's trinomials to each simplex coordinates ====%
%========== And Integrating all coefficients for this simplex ============%
h = waitbar(0,'Please wait');
for f = 1:nf
    waitbar(f/nf)

    j1 = face(f,2);
    j2 = face(f,3);
    j3 = face(f,4);
    
    x = [vertex(j1,2),vertex(j2,2),vertex(j3,2)];
    y = [vertex(j1,3),vertex(j2,3),vertex(j3,3)];
    z = [vertex(j1,4),vertex(j2,4),vertex(j3,4)];
    
    detJ = det([x;y;z]);
    
    Xn = Pascal(x,n);
    Yn = Pascal(y,n);
    Zn = Pascal(z,n);
    
    for nn = 0:n
        if nn == 5
            disp('ali')
        end
        for mm = 0:nn
            
            [cpri,spri] = Cprime(C(nn+1,mm+1).trinom,S(nn+1,mm+1).trinom,Xn,Yn,Zn);
         
            for i = 0:nn
                for j = 0:nn-i
                    
                    k = factorial(i)*factorial(j)*factorial(nn-i-j)/factorial(nn+3)/R^nn;
                    c(nn+1,mm+1) = c(nn+1,mm+1) + k*detJ*cpri(i+1,j+1);
                    s(nn+1,mm+1) = s(nn+1,mm+1) + k*detJ*spri(i+1,j+1);
                    
                end
            end
        end
    end
end
close(h)
c = c/volume;
s = s/volume;


%======================== Printimng the results ==========================%
fprintf('Volume :  %f(km^3) \n\n',volume)

for i = 0:n
    for j = 0:i
        fprintf('%d  %d :    %f     %f \n',i,j,c(i+1,j+1),s(i+1,j+1))
    end
end


%========================== Printimng the body ===========================%
plot3(vertex(:,2),vertex(:,3),vertex(:,4),'.','markersize',1)
grid on
