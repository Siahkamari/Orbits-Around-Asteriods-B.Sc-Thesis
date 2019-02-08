%================================ CPrime =================================%

% This function transforms the trinomials for one pair of the spherical
% harmonics coefficient to the coordinate frame of one simplex given the
% trinomial for those pair of coefficients, and 2 dimentional pascal
% triangul for (x1 x2 x3) & (y1 y2 y3) & (z1 z2 z3) which 1 , 2 , 3 each
% stands for one vertex of the simplex

function [cpri,spri] = Cprime(C,S,Xn,Yn,Zn)

n = numel(C(1,:))-1;
cpri = zeros(n+1,n+1);
spri = zeros(n+1,n+1);

for i = 0:n
    for j = 0:n-i
        Xi = Xn(i+1).trinom;
        Yj = Yn(j+1).trinom;
        Zk = Zn(n-j-i+1).trinom;
        cpri = cpri + C(i+1,j+1)*(Multiply(Multiply(Xi,Yj),Zk));
        spri = spri + S(i+1,j+1)*(Multiply(Multiply(Xi,Yj),Zk));
    end
end
        
end