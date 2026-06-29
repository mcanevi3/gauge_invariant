clear;clc;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
projR=@(C,R)eye(size(C,1))-C*pinv(C'*R*C)*C'*R;
C=[1,0;
   1,0;
   1,0;
   ];
% C=[1,0;
%    1,0;
%    1,0;
%    ];

syms f1 f2 f3 f4 f5 f6 real;
syms n1 n2 n3 real;
syms x1 x2 real;
x=[x1;x2];

y=C*x+[f1;f2;f3];



R=diag([10,1,100]);
Pival=projR(C,R);
vpa(Pival*y,5)
%[A,b]=equationsToMatrix(Pival*y)

Pival=proj(C);
vpa(Pival*y,5)
