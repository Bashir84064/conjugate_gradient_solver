%%%
%creates a 4096 x 4096 matrix for a 2d laplacian grid .
%applies CGM for three types of right hand sides i.e. b vector.
%1)b is uniformally distributed random vector.
%2)b is normally distributed random vector.
%3)b is a eigenvector corresponding to smallest eigen value.

%creation of matrix for 2d laplace operator for 64 unknowns in each operation 
internalPoints=64;
e   = ones(internalPoints,1);
spe = spdiags([e -2*e e], -1:1,...
   internalPoints,internalPoints);
Iz  = speye(internalPoints);
sp  = kron(Iz,spe)+kron(spe,Iz);
A=full(sp);

%4096x1  uniformally distributed random vector

xs = rand(4096,1);
b = A*xs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 10^-6;
m = 100;
%(1) Solving for b- uniformally distributed vector,and plotting the error and residue 
X= cg(A,b,tol,m);
error_matrix = (X - xs);
%norm of difference
error_norm = vecnorm(error_matrix);
%residue
r = zeros(length(xs),m);
r = b - A*X;
res_norm = vecnorm(r);
%axe  = 1:m;
plot(error_norm,'xr');
hold on ;
xlim([0,120]);
xlabel("number of iterations");
ylabel("error and residue");
plot(res_norm,'g');
legend("error 2-norm","residue 2-norm");
title("Application of CG on a 2D Laplace Operator with uniformally distributed random vector");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(2) Solving for b1-normally distributed vector,finding exact solution by
%matlab \ routine and then plotting the norm(x_cgs - A\b)
b1 = randn(4096,1);
X1 = cg(A,b1,tol,m);
x1 = A\b;
error1 = (X1 - x1);
error1_norm = vecnorm(error1);
figure();
plot(error1_norm,"ob");
hold on;

%(3) Solving for b2- eigen vector with smallest eigen value
[b2,small_eig] = eigs(A,1,'SM'); %GETTING THE eigenvector corresponding to smallest eigenvalue
X2 = cg(A,b2,tol,m);
error2 = (X2-x1);
error2_norm = vecnorm(error2);
plot(error2_norm,"xg");
xlabel("number of iterations");
ylabel("error ");

legend("normally distributed","eigenvector of smallest eigen value");
title("CGS with different b's");



