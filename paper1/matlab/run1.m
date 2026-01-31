clear;clc;
% States x v theta omega
B=1;
M=1;
Br=1;
J=1;
A=[[0,1;0,-B/M],zeros(2,2);zeros(2,2),[0,1;0,-Br/J]];
Bu=[[0;1/M],zeros(2,1);zeros(2,1),[0;1/J]];
Bw=zeros(5,1);
Dw=[1;1;1;1;1];
C=[1,0,0,0;0,1,0,0;0,1,0,0;0,0,1,0;0,0,0,1];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
proj_pi=proj(C);

proj_pi
proj_pi*C
proj_pi*null(C')
