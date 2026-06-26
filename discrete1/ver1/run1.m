clear;clc;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C'; 
projR=@(C,R)eye(size(C,1))-C*pinv(C'*R*C)*C'*R; 

A=[1,2;-3,-4];
Bu=[1;0];
Bw=[0.5;0.5];
C=[1,0];

ts=0.01;
Ad=eye(2)+ts*A;
Bud=ts*Bu;
Bwd=ts*Bw;
Cd=C;

O1=[Cd;Cd*Ad;Cd*Ad*Ad];
Huk=[0*Cd*Bud;Cd*Bud;Cd*Ad*Bud];
Fuk=[0*Cd;Cd;Cd*Ad];
Huk1=[0*Cd*Bud;0*Cd*Bud;Cd*Bud];
Fuk1=[0*Cd;0*Cd;Cd];

Pival=projR(O1,diag([1,1,1]));

syms x1 x2 real;
syms uk uk1 real;
syms fk_1 fk_2 fk1_1 fk1_2 real;

Yk=O1*[x1;x2]+Huk*uk+Fuk*[fk_1;fk_2]+Huk1*uk1+Fuk1*[fk1_1;fk1_2];
temp=Pival*(Yk-Huk*uk-Huk1*uk1);
for i=1:3
    disp(vpa(temp(i),4))
end
Pival*Fuk
Pival*Fuk1