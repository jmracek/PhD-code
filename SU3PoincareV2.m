function [ p ] = SU3PoincareV2( alpha )
%This program computes the equivariant Poincare polynomial of the loop
%group of SU(3) with the TxS1 action 

%How big do you want the simulation
N = 10;

%Initialize variables to use later
T = 0;
correct = [0 0 0];
syms t;

%Generators for root lattice
X1 = [1 0 0; 0 -1 0; 0 0 0];
X2 = [0 0 0; 0 1 0; 0 0 -1];

%Orthonormal basis for Cartan, K(x,y) = 0.5*tr(xy) Killing form
E1 = [1 0 0; 0 -1 0; 0 0 0];
E2 = [1 0 0; 0 1 0; 0 0 -2]/sqrt(3);

%Fundamental weights
L1 = [2/3 -1/3 -1/3];
L2 = [-1/3 2/3 -1/3];
L3 = [-1/3 -1/3 2/3];

%Iterate over Weyl invariant finite sets; these are indexed by the integral
%dominant weights

for m = 0:10
    for n = 0:10
        
        %Choice of integral dominant weight
        D = m*L1 - n*L3;
        
        %Compute the Weyl conjugates of D
        if (n == 0)||(m==0)
            %Adjoint orbit is a P2
            Delta = WeylOrbit(D,0);
        else
        end
    
        %Mod out the root lattice to determine proper correction
        r = mod(n-m,3);
        if r==0
            %You're already on the lattice. Don't do anything.
            correct = diag([0 0 0]);
        elseif r==1
            %You need to shift by a positive fundamental weight to get on
            %the root lattice
            correct = diag(L1);
        elseif r==2
            %You need to shift by a negative fundamental weight to get on
            %the root lattice
            correct = diag(-L1);
        end
               
        %Compute all the translates by roots
        for j = -N:N
            for k = -N:N
                
                %Translation amount
                T = j*E1+k*E2+correct;
                
                %Compute the translate of the Weyl orbit
                for ind = 1:length(Delta)
                    Delta{ind} = Delta{ind} + T;
                end
                
                %Figure out which \beta this orbit came from
                
                %Pick three points in the plane defining the image of the
                %critical set
                p1 = [ 0.5*trace(Delta{1}*E1) 0.5*trace(Delta{1}*E2) 0.5*trace(Delta{1}*Delta{1}) ];
                p2 = [ 0.5*trace(Delta{2}*E1) 0.5*trace(Delta{2}*E2) 0.5*trace(Delta{2}*Delta{2}) ];
                p3 = [ 0.5*trace(Delta{3}*E1) 0.5*trace(Delta{3}*E2) 0.5*trace(Delta{3}*Delta{3}) ];
                
                %The normal and shift of the plane are
                n = cross((p1-p2),(p1-p3));                
                D = dot(n,p1);
                
                gamma = @(t) alpha+n*t;
                
                %This is the \beta (closest point to the origin of a convex
                %subset of the weights)
                beta = solve(dot(gamma(t),n)-D==0);
                
                %If the z-component of beta is negative, then this
                %corresponds to a critical set of infinite index and we can
                %quit this iteration immediately
                
                Do this later
                
                %Now we need to solve the linear programming problem; is
                %the point beta on the interior of the convex hull of the
                %weights?
                
                for ind = 1:length(Delta)
                    
                end
                
                
                
            end
        end
        
        
        
    end
end



end