clear;clc;

% 2 DGUs Scenario 1
VDC=100;
Ctval=2.2*1e-3;
Ltval=1.8*1e-3;
Rtval=0.2;
fsw=10*1e3;
L12=1.8*1e-6;
R12=0.05;

C1=2.2*1e-3;
R1=0.2;
L1=1.8*1e-3;
C2=2.2*1e-3;
R2=0.2;
L2=1.8*1e-3;

% xdoti = Aii xi+ Bi ui+ Mi di+ zeta i
% yi=Ci xi
% zi=Hi yi
% xi=[Vi,Iti]^T, uti=Vti, di=ILi , zi=Vi
% zetai=Aij xj  coupling with DGU j

A=[[0,inv(C1);-inv(L1),-R1*inv(L1)], 
   [inv(R12*C1),0;0,0];
   [inv(R12*C2),0;0,0],
   [0,inv(C2);-inv(L2),-R2*inv(L2)] 
   ];
B=[
   [0;inv(L1)],   0*[0;inv(L1)];
   0*[0;inv(L2)],   [0;inv(L2)];
   ];
Bd=[
   [-inv(C1);0] ,  0*[-inv(C1);0];
   0*[-inv(C2);0],   [-inv(C2);0]
   ];

   
C=[1,0,0,0;
   1,0,0,0;
   0,1,0,0;
   0,1,0,0;
   0,0,1,0;
   0,0,0,1
];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);
