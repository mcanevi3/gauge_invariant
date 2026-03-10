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

syms l1 l2 l3 l4 l5 l6 real;
L=[l1,l2,l3;l4,l5,l6];

syms e1 e2 n1 n2 n3 real;
e=[e1;e2];
n=[n1;n2;n3];
edot=(A-L*C)*e+L*Dn*n
