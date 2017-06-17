function [ daplot,data1,data2 ] = plotorbits( funHandle )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%Pick how many one parameter subgroups to probe
m = 1;

%Choose number of steps;
NumSteps = 7;

%Generate the unit vectors in the Lie algebra for the 1-PS
vects = zeros(m,2);
for j=1:m
%     vects(j,:) = [real(exp(pi*1i*j./m/2)),imag(exp(pi*1i*j./m/2))];
vects(j,:) = [0 1];
end
vects

loopEnergy = zeros(2*NumSteps+1,m);
loopMomentum = zeros(2*NumSteps+1,m);

%For each one parameter subgroup
for j = 1:m
    
    tempLoop = funHandle;
    tempLoopN = funHandle;
    
    [loopEnergy(1,j),loopMomentum(1,j)] = EnergyMoment(funHandle);
    
    %Act the given 1PS on the loop and compute the momentum and
    %energy
    for k = 1:NumSteps
        tempLoop = act1PS(tempLoop,vects(j,:));
        tempLoopN = act1PS(tempLoopN,-vects(j,:));
        [loopEnergy(2*k,j),loopMomentum(2*k,j)] = EnergyMoment(tempLoop)
        [loopEnergy(2*k+1,j),loopMomentum(2*k+1,j)] = EnergyMoment(tempLoopN)
    end
end

ROYGBV = {'ro-','ko-','yo-','go-','bo-','mo-'};

x = [-7:1:7];
x2 = x.^2;

hold on
for j = 1:m
A = plot(real(loopMomentum(:,j)),real(loopEnergy(:,j)),ROYGBV{j});
end
A = plot(x,x2,'ko-')
hold off

daplot = A;
data1 = loopEnergy;
data2 = loopMomentum;

end

