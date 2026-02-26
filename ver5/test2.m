clear;clc;
%%% LO model1

A=-1;
Bw=1;

C=[1;2];
D=[3;3];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C)

lambda=-4;
l2=-lambda-4/3;
l1=1/3-l2;
L=[l1,l2];
eig(A-L*C)
Bw-L*D

