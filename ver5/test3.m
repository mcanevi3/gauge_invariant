clear;clc;

A=[1,2;
  -3,-4];
Bw=[1;
    1];
C=[2,3;
   1,2;
   3,2];
Dw=[1;
    1;
    -2];
Dn=eye(3);

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);

syms x1 x2 w n1 n2 n3;
y=C*[x1;x2]+Dw*w+[n1;n2;n3];

Pivalw=proj(Dw);

pinv(Pival*Dw)*Pival*y
Dw*w+[n1;n2;n3]

pinv(Dw)*(Dw*w+[n1;n2;n3])
