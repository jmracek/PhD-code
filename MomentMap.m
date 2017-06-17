function [p,E] = MomentMap(alpha,beta,t)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

EInt = 0;
PInt1 = 0;
PInt2 = 0;

E = 0;
p = 0;

ap = diff(alpha)./diff(t);
ap(length(t)) = ap(1);
bp = diff(beta)./diff(t);
bp(length(t)) = bp(1);

ENIG =  (abs(ap).^2 + abs(bp).^2)./(2*pi);
PIG1 =  (conj(alpha).*ap+conj(beta).*bp)./(2*pi);
PIG2 = (-ap.*beta+alpha.*bp)./(2*pi);

%Integrate using trapezoid method
for l = 1:length(t)-1
    EInt = EInt + (ENIG(l+1)-ENIG(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*ENIG(l);
    PInt1 = PInt1 + (PIG1(l+1)-PIG1(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*PIG1(l);
    PInt2 = PInt2 + (PIG2(l+1)-PIG2(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*PIG2(l);
end

E = EInt;
p = -0.5*trace([1i 0; 0 -1i]*[PInt1 -conj(PInt2); PInt2 conj(PInt1)]);

return

