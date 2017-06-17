function [B] = FourierLoop(ap,bp)
%This function is used to compute the Fouri
N = 6;
n = 1;
m = 1;

A = zeros(2*(2*N+1));

for j = 1:4:4*N
    k = -N+j-1;
    %order: a-k, b-k, ak, bk
    
    A(j,j) = bp;
    A(j,j+3) = 1i*2*n/m - 1i*k+ap;
    
    A(j+1,j+1) = 1i*2*n/m + 1i*k+ap;
    A(j+1,j+2) = bp;
    
    A(j+2,j+2) = ap-1i*k;
    A(j+2,j+1) = -bp;
    
    A(j+3,j) = ap+1i*k;
    A(j+3,j+3) = -bp;
end
A(length(A)-1,length(A)-1) = ap;
A(length(A)-1,length(A)) = -bp;
A(length(A),length(A)-1) = bp;
A(length(A),length(A)) = 1i*2*n/m+ap;

b = zeros(length(A),1);
coefs = A\b

t = [0:0.001:2*pi];

B = A;
end
