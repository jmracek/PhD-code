function [ Loop ] = AlmostComplex_v3(X,gamma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Number of Fourier coefficients
N = 40;
xi = zeros(2*N+1,4);

t = [0:2*pi/(length(X)-1):2*pi];

Y = X;
IG = zeros(length(X),1);
iter = 1;
temp = cell(length(X),1);

for i = 1:length(gamma)
   Y{i} = gamma{i}\X{i};
end

%Find the Fourier coefficients and apply the almost complex structure
for k = -N:N
    ind = k+N+1;
    %Finding the coefs
    for j1 = 1:2
        for j2 = 1:2
            
            Int = 0;
            
            %Initialize the integrand of the Fourier transform
            for l = 1:length(X)
                IG(l) = Y{l}(j1,j2).*exp(-1i.*k.*t(l))/(2*pi);
            end
            
            %Evaluate the integral using the trapezoid method
            for l = 1:length(t)-1
                Int = Int + (IG(l+1)-IG(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*IG(l);
            end
            
            %Store the Fourier coefficient
            xi(ind,iter) = Int;
            iter = iter+1;
            
        end
    end
    
    iter = 1;
    

    %Applying the ACS
    if (k > 0)
        xi(ind,:) = 1i*xi(ind,:);
    elseif (k < 0)
        xi(ind,:) = -1i*xi(ind,:);
    end
        
end

%Determine the zero'th Fourier mode
xi(N+1,:) = 0;
xi(N+1,:) = -sum(xi);


%Initialize a temporary variable for reconstructing the Fourier series
tempLoop = zeros(2,2);

%Rebuild the loop in Omega\frak{g}
for j = 1:length(X)
    for k = -N:N
        ind = k+N+1;
        tempLoop = tempLoop + [xi(ind,1) xi(ind,2); xi(ind,3) xi(ind,4)].*exp(1i.*k.*t(j));
    end
    
    temp{j} = tempLoop;
    tempLoop = zeros(2,2);
end


%Translate back to the tangent space above gamma
for j = 1:length(gamma)
   temp{j} = gamma{j}*temp{j};
end

Loop = temp;
end

