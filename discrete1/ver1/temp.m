clear;clc;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
C=[1;2];

Pival=eye(2)-C*pinv(C'*C)*C';
Pival

syms r1 r2 real;
Pival=eye(2)-C*pinv(C'*[r1,0;0,r2]*C)*C'*[r1,0;0,r2];

vpa(Pival,4)

R=[1,2;2,2];
Pival=eye(2)-C*pinv(C'*R*C)*C'*R;
Pival
