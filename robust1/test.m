clear;clc;

A=-1;
Bw=1;
Dw=[1;2]*1+[2;-1]*1;

C0=[1;2];
Ctheta=[2;3];

H=diag(C0);
syms d1 d2 real;
E=1;
H*[d1;d2]*E

theta_bar=1;


theta_test=-0.5;
C=C0+(theta_test)*Ctheta;
 
r=size(Bw,2);
n=size(A,1);
p=size(C,1);

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
DeltaM = Y*Ctheta;
[ A'*P+P*A - Y*C0 - (Y*C0)' + theta_bar^2*S, ...
  P*Bw - Y*Dw, Cz', DeltaM;
  (P*Bw - Y*Dw)', -gmma*eye(r), zeros(r,n), zeros(r,n);
  Cz, zeros(n,r), -gmma*eye(n), zeros(n,n);
  DeltaM', zeros(n,r), zeros(n,n), -S ] <= 0;
cvx_end;

L=inv(P)*Y;
L
L=[10,10];
eig(A-L*C)

Pi0=proj(C0);
Piub=proj(C0+theta_bar*Ctheta);
Pitheta=proj(Ctheta);

syms alpha real;
C=C0+(alpha)*Ctheta;

xt=0.01;
yt=C*xt;

temp=Pi0*yt;
pinv(Pi0*Ctheta)*temp
