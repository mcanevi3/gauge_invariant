clear;clc;

A=-1;
Bw=1;

C=[1;2];
D=[3;3];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
proj(C)
