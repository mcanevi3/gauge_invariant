clear;clc;

% 2 DGUs Scenario 1
VDC=100;
Ctval=2.2*1e-3;
Ltval=1.8*1e-3;
Rtval=0.2;
fsw=10*1e3;
L12=1.8*1e-6;
R12=0.05;

Ct=[Ctval,Ctval];
Lt=[Ltval,Ltval];
Rt=[Rtval,Rtval];
% xdoti = Aii xi+ Bi ui+ Mi di+ zeta i
% yi=Ci xi
% zi=Hi yi
% xi=[Vi,Iti]^T, uti=Vti, di=ILi , zi=Vi
% zetai=Aij xj  coupling with DGU j

A11=@(i)[-inv(R12*Ct(i)),inv(Ct(i));-inv(Lt(i)),-Rt(i)*inv(Lt(i))];
A12=@(i)[inv(R12*Ct(i)),0;0,0];
B1=@(i)[0;inv(Lt(i))];
M1=@(i)[-inv(Ct(i));0];
C1=eye(2);
H1=[1,0];

A=[A11(1),A12(1);A12(2),A11(2)];
B=[B1(1),0*B1(1);B1(2)*0,B1(2)];
M=[M1(1),0*M1(1);M1(2)*0,M1(2)];

C=[1,0,0,0;
   1,0,0,0;
   0,1,0,0;
   0,1,0,0;
   0,0,1,0;
   0,0,0,1
];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);
