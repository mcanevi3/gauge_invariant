clear;clc;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';

C=[1,2;1,2;1,2;1,2;1,2];
Dw=[1;1;1;1;1];

disp(string(size([C,Dw],1))+" eqns");
disp("rank:"+string(rank([C,Dw])));
syms x1 x2 real;
syms w real;
syms f1 f2 f3 f4 f5 f6 real;

%y=C*[x1;x2]+Dw*w+[f1;f2;f3;f4;f5;f6];

Pi1=proj(C);
Pi2=proj(Pi1*Dw);

T=Pi2*Pi1;
temp=T*[f1;f2;0;0;0];
vpa(temp,3)

temp=T*[rand();rand();0;0;0];
temp
%[U,S,V]=svd(T);
%S

