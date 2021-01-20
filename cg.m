
%CGS algorithm 
%Inputs : matrix A ,vector b,tolerance and number of iterations
%Ouput :  matrix X,solution of Ax = b where every column of X is an iterate solution
%stops :until convergence condition is met or number of iterations are over

function X = cg(A,b,tol,m)
n= size(A,1);
X = zeros(n,m); %already caters for x_0 = [0,0,0....,0]
r = zeros(n,m);
p = zeros(n,m);
%defining r0,p0
r(:,1) = b - A*X(:,1);
p(:,1) = r(:,1);
%norm r0
norm_r0 = norm(r(:,1));
%breaking condition
if norm_r0 < tol
     return;
end

%for loop
for j = 1:m
    %only one matrix vector product in every iteration
    mat_vect = A * p(:,j);
    gama = (r(:,j)' * r(:,j))/(p(:,j)' * mat_vect);
    X(:,j+1) = X(:,j) + gama * p(:,j);
    r(:,j+1) = r(:,j) - gama * mat_vect;
    %breaking condition
    if (norm(r(:,j+1))/norm_r0)< tol
        %fprintf("CGs Converged at %d steps /n",j);
        break;
    end
    
    conj_direction = (r(:,j+1)' * r(:,j+1))/(r(:,j)' * r(:,j));
    p(:,j+1) = r(:,j+1) + conj_direction * p(:,j);
    
end




end
