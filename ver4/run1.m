clear;clc;

% A=-1;
% Bw=[1];
% C=[2;1];

A=[1,2;-3,-4];
Bw=[1;1];
C=[2,1;1,2];
n=2;
p=2;

cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(p,n);
variable ro;
minimize(ro);
P>=0;
[A'*P+P*A-Y*C-C'*Y',P*Bw-Y*C;(P*Bw-Y*C)',-ro*eye(1)]<=0;
cvx_end;

L=inv(P)*Y;
gmma=sqrt(ro);

disp(eig(A-L*C))
% G=ss(A-L*C,(Bw-L)*C,1,0);
% tf(G)