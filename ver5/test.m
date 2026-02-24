clear;clc;

A=-1;
Bw=1;

C0=[1;2];
H=[1,0;0,2];
E=1;

D0=[3;3];
G=[1,0;0,1];
F=1;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';

syms d1 d2 d3 d4 real;
syms x w real;

%y=(C0+H*[d1;d2]*E)*x+(D0+G*[d3;d4]*F)*w;
%Pival0=proj(C0);

PivalC=proj(C0);
PivalD=proj(D0);
y=(C0+H*[d1;d2]*E)*x+D0*w;

PivalCD=proj(PivalD*C0)
PivalD
PivalCD*PivalD

