%============================ Coefficients ===============================%

% This function calculates the trinomials for spherical harmonics
% coefficients till degree n

function [C , S] = Coefficients(n)

C(n+1,n+1) = struct('trinom',[]);
S(n+1,n+1) = struct('trinom',[]);

x = [0  0
    1   0];
y = [0  1
    0   0];
z = [1  0
    0   0];
r2 = [1	0   1
    0	0   0
    1   0   0];

%=========================== Initializiation =============================%
C(0+1,0+1).trinom = 1;
S(0+1,0+1).trinom = 0;

if n ~=0
%     C(1+1,1+1).trinom = x/sqrt(3);
%     S(1+1,1+1).trinom = y/sqrt(3);
    
    % for Unnormalized coefficients
    C(1+1,1+1).trinom = x;
    S(1+1,1+1).trinom = y;
end


%============================== Sectorial ================================%
for i = 2:n
    
%     k = (2*i - 1)/sqrt(2*i*(2*i + 1));
    
    % for Unnormalized coefficients
    k = 1/(2*i);
    
    C(i+1,i+1).trinom = ...
        k*Multiply(x,C(i-1+1,i-1+1).trinom) -...
        k*Multiply(y,S(i-1+1,i-1+1).trinom);
    
    S(i+1,i+1).trinom = ...
        k*Multiply(y,C(i-1+1,i-1+1).trinom) +...
        k*Multiply(x,S(i-1+1,i-1+1).trinom);
    
end


%============================== Subdiagonal ==============================%
for i = 1:n
    
%     k = (2*i-1)/sqrt(2*i+1);
    
    % for Unnormalized coefficients
    k = 1;
    
    C(i+1,i-1+1).trinom =...
        k*Multiply(z,C(i-1+1,i-1+1).trinom);
    
    S(i+1,i-1+1).trinom =...
        k*Multiply(z,S(i-1+1,i-1+1).trinom);
    
end


%=============================== Vertical ================================%
for j = 0:n-2
    for i = j+2:n
        
%         k1 = (2*i-1)*sqrt((2*i-1)/((2*i+1)*(i+j)*(i-j)));
%         k2 = -sqrt((2*i-3)*(i+j-1)*(i-j-1)/((2*i+1)*(i+j)*(i-j)));
        
        % for Unnormalized coefficients
        k1 = (2*i-1)/(i+j);
        k2 = -(i-j-1)/(i+j);
        
        C(i+1,j+1).trinom = ...
            k1*Multiply(z,C(i-1+1,j+1).trinom) +...
            k2*Multiply(r2,C(i-2+1,j+1).trinom);
        
        S(i+1,j+1).trinom = ...
            k1*Multiply(z,S(i-1+1,j+1).trinom) +...
            k2*Multiply(r2,S(i-2+1,j+1).trinom);
        
    end
end

end

