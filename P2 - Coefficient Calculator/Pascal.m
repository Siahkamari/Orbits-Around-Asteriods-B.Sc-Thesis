%============================== Pascal ===================================%

% This function calculates the 2D Pascal triangul

function Y = Pascal(X,n)

Y(n+1) = struct('trinom',[]);
Y(1).trinom = 1;

for i = 1:n
    
    Y(i+1).trinom = zeros(i+1,i+1);
    Y(i+1).trinom(1:i,1:i) =  X(3)*Y(i).trinom;
    Y(i+1).trinom(2:i+1,1:i) = Y(i+1).trinom(2:i+1,1:i) + X(1)*Y(i).trinom;
    Y(i+1).trinom(1:i,2:i+1) = Y(i+1).trinom(1:i,2:i+1) + X(2)*Y(i).trinom;

end

end