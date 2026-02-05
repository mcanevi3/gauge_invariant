clear;clc;
rng(233)
% is there a solution Bw-L*Dw\neq 0 ??? Symbolic

n=2; % num of states
r=1; % size of w
p=3; % num of sensors
A=rand(n,n);
Bw=rand(n,r);
C=rand(p,n);

Cperp=null(C');
Pival=eye(p)-C*pinv(C'*C)*C';

Dw=rand(p,r);

L=sym('x',[n,p]);
expr=Bw-L*Dw;
disp("Bw-L*Dw:");
disp(expr);
sol=vpasolve(expr==0);
sol