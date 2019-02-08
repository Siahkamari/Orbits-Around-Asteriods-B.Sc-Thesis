%============================== Multiply =================================%

% This function multiplys two trinomials

function Y = Multiply(A,B)

nA = numel(A(:,1))-1;
nB = numel(B(:,1))-1;
n = nA + nB ;
Y = zeros(n+1,n+1);

for i = 0:nB
    for j = 0:nB - i
        
    Y(i+1:nA+i+1,j+1:nA+j+1) =  Y(i+1:nA+i+1,j+1:nA+j+1) + B(i+1,j+1)*A;

    end
end

end