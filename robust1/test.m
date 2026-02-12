clear;clc;

A=-1;
Bw=1;
Dw=[1;2]*1+[2;-1]*1;

C0=[1;2];
H=diag(C0);
E=1;

theta_bar=1;
 
r=size(Bw,2);
n=size(A,1);
p=size(C0,1);

proj=@(C)eye(p)-C*pinv(C'*C)*C';


Cz=eye(n);
cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable S(n,n) symmetric
variable Y(n,p); % n,p
variable gmma;
minimize(gmma);
P>=1e-5*eye(n);
S >= 1e-5*eye(n);
DeltaM = Y*C0;
[ A'*P+P*A - Y*C0 - (Y*C0)' + theta_bar^2*S, ...
  P*Bw - Y*Dw, Cz', DeltaM;
  (P*Bw - Y*Dw)', -gmma*eye(r), zeros(r,n), zeros(r,n);
  Cz, zeros(n,r), -gmma*eye(n), zeros(n,n);
  DeltaM', zeros(n,r), zeros(n,n), -S ] <= 0;
cvx_end;

L=inv(P)*Y;
L
L=[20,20];
eig(A-L*C0)

Pi0=proj(C0);

