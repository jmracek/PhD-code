function tube = GradFlow(gamma)
%This function inputs a loop \gamma \in \Omega SU(2) and computes its flow under the gradient of the norm squared of the moment map for the TxS1 action
%Here, \gamma is an array of 2x2 matrices storing the pointwise values of the loop
%The gradient flow equations are solved numerically using the implicit Euler method

%INPUT: gamma is a cell array of matrices representing a path in SU(2).  Each cell is the value of the path at a particular time interval between 0 and 2pi

shift = [-0.5 0.457];
NUM = 1;
tol = 1;

%cell array for storing final results
H = cell(1,length(gamma));
tempPath = cell(length(gamma),1);

for j = 1:length(gamma)
    H{1,j} = gamma{j};
    H{2,j} = gamma{j};
end

alpha = zeros(length(gamma),1);
beta = zeros(length(gamma),1);
EInt = 0;
PInt1 = 0;
PInt2 = 0;

Sz = [1i 0; 0 -1i];

Xt = cell(length(gamma),1);
Xs1 = cell(length(gamma),1);
Xmu = cell(length(gamma),1);

%step size
h = 2*pi/(length(gamma)-1);

t = [0:h:2*pi]';


for k=1:20
    k
    
    for k2 = 1:5
        k2
        % Calculate X_{\mu(\gamma)}
        % calculate E and p
        
        for j = 1:length(gamma)
            alpha(j) = H{k+1,j}(1,1);
            beta(j) = H{k+1,j}(2,1);
            tempPath{j} = H{k+1,j};
        end
        
        
        ap = diff(alpha)./diff(t);
        ap(length(t)) = ap(1);
        bp = diff(beta)./diff(t);
        bp(length(t)) = bp(1);
        
        ENIG =  (abs(ap).^2 + abs(bp).^2)./(4*pi);
        PIG1 =  (ap.*conj(alpha)+beta.*conj(bp))./(2*pi);
        PIG2 = (conj(alpha).*bp-beta.*conj(ap))./(2*pi);
        
        %Integrate using trapezoid method
        for l = 1:length(t)-1
            EInt = EInt + (ENIG(l+1)-ENIG(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*ENIG(l);
            PInt1 = PInt1 + (PIG1(l+1)-PIG1(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*PIG1(l);
            PInt2 = PInt2 + (PIG2(l+1)-PIG2(l))*(t(l+1)-t(l))/2 + (t(l+1)-t(l))*PIG2(l);
        end
        
        E = EInt;
        P = -0.5*trace([1i 0; 0 -1i]*[PInt1 -conj(PInt2); PInt2 conj(PInt1)]);
        
        EInt = 0;
        PInt1 = 0;
        PInt2 = 0;
        
        %Calculate the torus and S1 H.V.F along gamma
        for j=1:length(gamma)
            
            Xt{j} = Sz*H{k+1,j}-H{k+1,j}*Sz;
            Xs1{j} = [ap(j) -conj(bp(j)); bp(j) conj(ap(j))] - H{k+1,j}*[ap(1) -conj(bp(1)); bp(1) conj(ap(1))];
            
            Xmu{j} = Xt{j}.*(P-shift(1)) + Xs1{j}.*(E-shift(2));
            
        end
        
        %Solve the flow equations for each fixed t using backward Euler
        
        for j = 1:length(gamma)
            
            JV = AlmostComplex_v3(Xmu,tempPath);
            
            H{k+1,j} = H{k,j}+h.*JV{j};
        end
    end
    %Initial guess for the loop at the next time step
    for j = 1:length(gamma)
        H{k+2,j} = H{k+1,j};
    end
end

tube = H;

return

