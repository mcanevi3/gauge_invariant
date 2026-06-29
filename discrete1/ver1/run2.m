clear;clc;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C'; 

A=[1,2;-3,-4];
Bu=[1;0];
Bw=[0.5;0.5];
C=[1,0];
Dw=0.1;

ts=0.01;
Ad=eye(2)+ts*A;
Bud=ts*Bu;
Bwd=ts*Bw;
Cd=C;
Dwd=Dw;

Pival=proj([Cd;Cd*Ad;Cd*Ad*Ad]);

N=3;
x1=sym("x1",[1,N*2],"real");
x2=sym("x2",[1,N*2],"real");
w=sym("w",[1,N*2],"real");
n=sym("n",[1,N*2],"real");

y1=Cd*[x1(1:3);x2(1:3)]+Dwd*w(1:3)+n(1:3);
y2=Cd*[x1(2:4);x2(2:4)]+Dwd*w(2:4)+n(2:4);
y3=Cd*[x1(3:5);x2(3:5)]+Dwd*w(3:5)+n(3:5);

Yval=[y1;y2;y3];
Yval
vpa(Pival*Yval(:,1),5)