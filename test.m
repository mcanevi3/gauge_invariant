clear;clc;

A=rand(1,1);
Bw=rand(1,1);
C=rand(2,1);
Dw=[1;1];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);

Cperp=null(C');

Pival*C
Pival*Cperp

pinv(Pival*Dw)*Pival*Dw
