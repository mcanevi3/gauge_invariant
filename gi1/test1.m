clear;clc;

n=5;
C=ones(n,1);
lambda=eig(C*C');
lambda=lambda(lambda>0);
rho=max(lambda);
Pim=eye(n,n)-inv(rho)*C*C';

b=randn(n-1,1);
b=[b;-sum(b)];

disp("Pim*C:");
disp(Pim*C);
disp("b:");
disp(b);

disp("b'*C:");
disp(b'*C);
