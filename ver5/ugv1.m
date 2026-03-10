clear;clc;

B=1;
Br=1;
M=1;
J=1;
As=[[0,1;0,-B/M],zeros(2,2);zeros(2,2),[0,1;0,-Br/J]];
Bs=[0;1/M;0;1/J];
Cs=[1,0,0,0;0,1,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
Ds=eye(5);
%L=place(As',Cs',[-1,-2,-3,-4])';

Cz=[1,0,0,0];
n=size(As,1);
p=size(Cs,1);

cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,p); % n,p
variable gmma;
minimize(gmma);
P>=1e-3*eye(n);
[As'*P+P*As-Y*Cs-Cs'*Y',-Y*Ds,Cz';
(-Y*Ds)',-gmma*eye(p),zeros(p,1);
Cz,zeros(1,p),-gmma*eye(1)]<=0;
cvx_end;

L=inv(P)*Y;
disp(eig(As-L*Cs))
G=ss(As-L*Cs,(Bs-L)*Cs,eye(n),zeros(n,n));

syms b1 b2 b3 b4 b5 real;
L*Ds*[b1;b2;b3;b4;b5]
