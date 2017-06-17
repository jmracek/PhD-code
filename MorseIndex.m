function [ EV,index ] = MorseIndex(n,m,C)
%This function calculates the Morse index of critical points for the norm
%squared of the moment map of \Omega SU(2)

%Input: A path \gamma \in \Omega SU(2) which is fixed under some subtorus
%of the T \times S1 action

%Output: \nu(\gamma) the Morse index at gamma

%Number of positive Fourier modes to consider
N = 10;
%Initialize the Hessian
Hess = zeros(6*N);
%Initialize array for storing calculations
H = zeros(2);

time = 0:2*pi/50:2*pi;

%Initialize the path cell
gamma = cell(length(time),1);

% Initialize the deformation cell arrays
X = cell(length(gamma),1);
Y = cell(length(gamma),1);

h = 2*pi/(length(gamma)-1)/5;

tempPath1 = cell(length(gamma),1);
tempPath2 = cell(length(gamma),1);
tempPath3 = cell(length(gamma),1);
tempPath4 = cell(length(gamma),1);

%Calculate weights to join
k1 = n/m + C;
k2 = n/m - C;
%Let the path be halfway between the weights
A = -n/m;
B = C;
%This is the critical path
L = @(t) [(-k2-A)/(k1-k2)*exp(-1i*k1*t)+(A+k1)/(k1-k2)*exp(-1i*k2*t), -B/(k1-k2)*exp(-1i*k1*t)+B/(k1-k2)*exp(-1i*k2*t); B/(k1-k2)*exp(1i*k1*t)-B/(k1-k2)*exp(1i*k2*t), (-k2-A)/(k1-k2)*exp(1i*k1*t)+(A+k1)/(k1-k2)*exp(1i*k2*t)];
L2 = @(t) expm([1i 0; 0 -1i].*t);

for l = 1:length(gamma)
    gamma{l} = L2(time(l));
end

indj = 0;
indk = 0;

%Pick a pair of variations
for j = 0:6*N-1
    
    %Initizlize the deformation X
    if mod(j,6) == 0
        indj = indj+1
        for l = 1:length(time)
            X{l} =  h.*[1 0; 0 -1].*exp(1i.*indj.*time(l))+h*[-1 0; 0 1].*exp(-1i.*indj.*time(l));
        end
    elseif mod(j,6) == 1
        for l = 1:length(time)
            X{l} =  h*[0 1; 0 0].*exp(1i.*indj.*time(l))+h*[0 0; -1 0].*exp(-1i.*indj.*time(l))+h*[0 -1; 1 0];
        end
    elseif mod(j,6) == 2
        for l = 1:length(time)
            X{l} =  h*[0 0; 1 0].*exp(1i.*indj.*time(l))+h*[0 -1; 0 0].*exp(-1i.*indj.*time(l))+h*[0 1; -1 0];
        end
	elseif mod(j,6) == 3
        for l = 1:length(time)
            X{l} =  1i*h.*[1 0; 0 -1].*exp(1i.*indj.*time(l))-1i*h*[-1 0; 0 1].*exp(-1i.*indj.*time(l))-2i.*h.*[1 0; 0 -1];
        end
    elseif mod(j,6) == 4
        for l = 1:length(time)
            X{l} =  1i*h*[0 1; 0 0].*exp(1i.*indj.*time(l))-1i*h*[0 0; -1 0].*exp(-1i.*indj.*time(l))-2i.*h*[0 1; 1 0];
        end
    elseif mod(j,6) == 5
        for l = 1:length(time)
            X{l} =  1i*h*[0 0; 1 0].*exp(1i.*indj.*time(l))-1i*h*[0 -1; 0 0].*exp(-1i.*indj.*time(l))-2i*h*[0 1; 1 0];
        end
    end
    
    for k = 0:6*N-1
        
        %Initizlize the deformation Y
        if mod(k,6) == 0
            indk = indk + 1;
            for l = 1:length(time)
                Y{l} =  h*[1 0; 0 -1].*exp(1i.*indk.*time(l))+h*[-1 0; 0 1].*exp(-1i.*indk.*time(l));
            end
        elseif mod(k,6) == 1
            for l = 1:length(time)
                Y{l} =  h*[0 1; 0 0].*exp{1i.*indk.*time(l))+h*[0 0; -1 0].*exp(-1i.*indk.*time(l))+ h.*[0 -1; 1 0];
            end
        elseif mod(k,6) == 2
            for l = 1:length(time)
                Y{l} =  h*[0 0; 1 0].*exp(1i.*indk.*time(l))+h*[0 -1; 0 0].*exp(-1i.*indk.*time(l))+h.*[0 1; -1 0];
            end
elseif mod(k,6) == 3
        for l = 1:length(time)
            Y{l} =  1i*h.*[1 0; 0 -1].*exp(1i.*indk.*time(l))-1i*h*[-1 0; 0 1].*exp(-1i.*indk.*time(l))-2i.*h.*[1 0; 0 -1];
        end
    elseif mod(k,6) == 4
        for l = 1:length(time)
            Y{l} =  1i*h*[0 1; 0 0].*exp(1i.*indk.*time(l))-1i*h*[0 0; -1 0].*exp(-1i.*indk.*time(l))-2i.*h*[0 1; 1 0];
        end
    elseif mod(k,6) == 5
        for l = 1:length(time)
            Y{l} =  1i*h*[0 0; 1 0].*exp(1i.*indk.*time(l))-1i*h*[0 -1; 0 0].*exp(-1i.*indk.*time(l))-2i*h*[0 1; 1 0];
        end
    end
        end

        for iter=1:length(gamma)
            tempPath1{iter} = gamma{iter} - gamma{iter}*X{iter} - gamma{iter}*Y{iter};
            tempPath2{iter} = gamma{iter} - gamma{iter}*X{iter} + gamma{iter}*Y{iter};
            tempPath3{iter} = gamma{iter} + gamma{iter}*X{iter} - gamma{iter}*Y{iter};
            tempPath4{iter} = gamma{iter} + gamma{iter}*X{iter} + gamma{iter}*Y{iter};
        end
        
        H(1,1) = NormSquared(tempPath1);
        H(1,2) = NormSquared(tempPath2);
        H(2,1) = NormSquared(tempPath3);
        H(2,2) = NormSquared(tempPath4);
        
        Hess(j+1,k+1) = (H(1,1)-H(1,2) - H(2,1) + H(2,2))/(4*h^2);
        
    end
    indk = 0;
end
%

Hess
%Calculate the eigenvalues of the Hessian
EV = eig(Hess);

%Count the number of negative eigenvalues
countneg = 0;
countzero = 0;

for l = 1:length(EV)
    if EV(l) < -1e-3
        countneg = countneg + 1;
    elseif abs(EV(l)) < 1e-3
        countzero = countzero + 1;
    end
end

index(1) = countneg;
index(2) = countzero;

end
