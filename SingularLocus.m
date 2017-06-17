function [ A ] = SingularLocus()

%This function plots the set of singular values of the moment map \mu: \Omega SU(2) \to Lie(TxS1).  It is the collection of all straight lines joining
%the pairs of points of the form (n,n^2/2)
%
N = 30;
syms q;
P = gcf;
hold on
for j = -N:N
   for k = j+1:N
       L = @(q) (j+k)*q-j*k;
       fplot(L,[min([j k]),max([j k])])
   end
end
hold off


end
