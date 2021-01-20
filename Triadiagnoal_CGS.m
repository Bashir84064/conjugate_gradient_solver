% Plots the convergence of CGS of a tridiagnoal matrix of size 
% 1) 10x10
% 2) 100 x 100
% 3) 1000 x 1000
%plot number of iterations along the eucledian norm of exact solution and
%approximated solution.Algorithm reaches convergence at N/2 steps as shown
%in all three figures.



tol = 10^-6;
main = 2;
b = -1;
c = -1;
N =[10,100,1000];
name = ["10 iterations","100 iterations","1000 iterations"];
style = ["og","or","oy"];
for i = 1:length(N)

%A = diag(main*ones(1,N(i))) + diag(b*ones(1,N(i)-1),1) + diag(c*ones(1,N(i)-1),-1);
A = full(gallery('tridiag',N(i),-1,2,-1));
b=zeros(N(i),1);
b(1) = 1;
b(end) = 1;

x_exact = ones(N(i),1);
X =  cg(A,b,tol,N(i));
% difference between exact and approximates solution iterates
% x(1),x(2),,,x(j)
e_matrix = (X - x_exact);
%norm of difference
euc_norm = vecnorm(e_matrix);



%PLOTTING
figure(i);    
axe  = 1:N(i);
plot(axe,euc_norm,style(i));
title(name(i));


xlabel("No of iterations ");
 ylabel(" eucledian norm of solution");

end
