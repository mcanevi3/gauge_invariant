clear;clc;

A=-1;
Bw=1;
Dw=[1;2]*1+[2;-1]*1;

C0=[1;2];

r=size(Bw,2);
n=size(A,1);
p=size(C0,1);


proj=@(C)eye(p)-C*pinv(C'*C)*C';
Pi0=proj(C0);

syms d1 d2 real;
H=diag(C0);
E=1;

Ctheta=H*[d1;d2]*E;
Pi0*Ctheta 

Ctheta=H*[1;1]*E;
Pi0*Ctheta 

Ctheta=H*[1;1+1e-2]*E;
Pi0*Ctheta 

% if H=I and d1=d2 then Pi*(C0+H[d1;d2]E)=0
% otherwise recoverable

