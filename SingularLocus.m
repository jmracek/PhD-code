function [ A ] = SingularLocus()
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
