clear;clc;

A=[1,2;-3,-4];
B=[1;0];
C=[1,0];
n=size(A,1);
p = size(C, 1);

rho=inv(1000)^2;
cvx_begin sdp quiet;
variable P(n,n) symmetric
variable Gamma(p,p) symmetric
variable Y(n, p)
variable rhon
minimize(rhon)
Gamma>=1e-5;
P >= 1e-6 * eye(n); 
[A'*P+P*A-Y*C-C'*Y'+eye(n),C'*rho*Gamma'-Y,-Y;
 rho*Gamma*C-Y',rho*Gamma+rho*Gamma',rho*Gamma;
 -Y',rho*Gamma',-rhon*eye(1)]<=0;
cvx_end 

if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    L = inv(P)*Y;

    disp("L");
    disp(L);
    disp("eig(A-LC):");
    disp(eig(A-L*C));
    disp("Gamma");
    disp(Gamma);
else 
    disp("Status:"+string(cvx_status));
end