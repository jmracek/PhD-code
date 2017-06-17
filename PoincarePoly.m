function [ A ] = PoincarePoly(alpha)
%This program computes the Poincare polynomial of the symplectic quotient
%of \Omega SU(2) by the TxS1 action, reduced at a (generically chosen)
%level alpha.
syms p;
%Number of 1PS to include in the calculation
N = 50;
hold on
plot(alpha(1),alpha(2),'r*')
hold off
%Initialize the Poincare polynomial
Pq = @(q) 0;


%Find contributions to the Poincare polynomial coming from subtori as
%stabilizers
for n = -N:N

    if n == 0
        amin = ceil(abs(alpha(1)));
    for a = amin:amin+N
index = TbIndex(a,[0,1]);
        if index ~= pi
        Pq = @(q) Pq(q) + q^(index)./(1-q^2);
        end
    end

    else
    %Fix a slope, corresponding to a particular class of subtori
    MoN = 2./n;
    L = @(p) MoN.*p+alpha(2)-MoN.*alpha(1);
    %Find where the line intersects the parabola E = p^2;
    SOLS = solve(L(p) == p.^2);
    SOLS = double(SOLS);
    
    if imag(SOLS(1))~=0
        %If the line doesn't intersect the parabola, then this 1PS doesn't
        %contribute anything.
        continue;
    else
        amin = ceil(min(SOLS));
        amax = floor(max(SOLS));
        
    end
    
    
    
    %Find where the image of the fixed point sets of T_\beta intersects the
    %line L_beta
    for a = amin:amax
        %Figure out which two weights are on im(T_\beta)
        C = sqrt((-a-n./2).^2);

        a1 = -n./2+C;
        a2 = -n./2-C;
        
        L2 = @(p) (a1+a2)*p-a1*a2;
        PT = solve(L(p)==L2(p));
        PT = double(PT);
        y = L2(PT);
        %Determine whether the solution PT is inside the convex hull
        B = floor(PT);
        T = ceil(PT);
        Ltemp = @(p) p.*(T.^2-B.^2)./(T-B)-T.*B;
        ytemp = Ltemp(PT);
        if ( ytemp > y)
            %Then the point PT is not inside the convex hull
            continue;
        else
            %Calculate the index and add a contribution to Pq 
            
            %This part was the distance from the origin (alph1,alph2) to
            %the point of intersection.
            beta(1) = PT-alpha(1);
            beta(2) = y-alpha(2);
            
            if ((beta(1) >= 0)||(beta(1) <= 0))&&(beta(2) >= 0)
                index = TbIndex(a1,[n,2]);
                if index ~= pi
                    Pq = @(q) Pq(q) + q^(index)./(1-q.^2);
                end
%                 if index == 0
%                     n
%                 end
            end

 
        end
        
    end
end
    
end

%Find contributions to Pq from the fixed points
for l = -N:N
    index = TbIndex(l,[l-alpha(1),l.^2-alpha(2)]);
    if index == pi
        continue;
    else
        Pq = @(q) Pq(q) + q^(index)./(1-q.^2).^2;
    end
    if index ==0
        l
    end
end

A = Pq;
end

