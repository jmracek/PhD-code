function [ mu2 ] = NormSquared( gamma )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

t = [0:2*pi/(length(gamma)-1):2*pi];
E = 0;
p = 0;
EInt = 0;
PInt1 = 0;
PInt2 = 0;

eps = 0.99;

for j = 1:length(gamma)
    alpha(j) = gamma{j}(1,1);
    beta(j) = gamma{j}(2,1);
end

ap = diff(alpha)./diff(t);
ap(length(t)) = ap(1);
bp = diff(beta)./diff(t);
bp(length(t)) = bp(1);

ENIG =  (abs(ap).^2 + abs(bp).^2)./(2*pi);
PIG1 =  (ap.*conj(alpha)+beta.*conj(bp))./(2*pi);
PIG2 = (conj(alpha).*bp-beta.*conj(ap))./(2*pi);

%Integrate using trapezoid method
for l = 1:length(t)-1
    EInt = EInt + (ENIG(l+1)-ENIG(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*ENIG(l);
    PInt1 = PInt1 + (PIG1(l+1)-PIG1(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*PIG1(l);
    PInt2 = PInt2 + (PIG2(l+1)-PIG2(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*PIG2(l);
end

E = EInt;
p = -0.5*trace([1i 0; 0 -1i]*[PInt1 -conj(PInt2); PInt2 conj(PInt1)]);

mu2 = p;

end

