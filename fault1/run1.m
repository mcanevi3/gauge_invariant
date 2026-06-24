clear;clc;

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
C=[2,1;
   2,1;
   2,1;
   2,1;
   ];

Pival=proj(C);

syms f1 f2 f3 f4 f5 f6 real;
syms x1 x2 real;
x=[x1;x2];

%y=C*x+eye(6)*[0;0;0;1;0;0];
%vpa(Pival*y,5)

n=rand(4,1);
n=n/sqrt(var(n));
n=n-mean(n);
n=1*n;

y=C*x+[f1;f2;f3;f4]+n;
vpa(Pival*y,5)

y=C*x+[f1;0;0;0]+n;
vpa(Pival*y,5)

y=C*x+[0;f2;0;0]+n;
vpa(Pival*y,5)

y=C*x+[0;0;f3;0]+n;
vpa(Pival*y,5)

y=C*x+[0;0;0;f4]+n;
vpa(Pival*y,5)

% proj(kron(randn(1,2),ones(5,1)))*[f1,0;0,f2;0,0;0,0;0,0]

Pival*[0;rand();0;0]