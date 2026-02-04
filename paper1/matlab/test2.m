clear;clc;
%%% solvibg for L and solving for Lmbda
A=-1;
Bw=1;
C=[1;2];
Dw=[1;2];
% L=[1,1];
n=size(A,1);
p=size(C,1);

angle=acosd(Dw'*C/(norm(C)*norm(Dw)));
fprintf("Angle: %f\n",angle);

% A=[1,2;-3,-4];
% Bw=[1;1];
% C=[1,2;2,4;1,3];
% Dw=[1;1;1];
% L=place(A',C',[-2,-3])';
% eig(A-L*C)

Pival=eye(p)-C*pinv(C'*C)*C';

% cvx_clear;
% cvx_begin sdp quiet;
% variable Lmbda(n,p)
% minimize(norm(Bw+(-L-L*Pival+Lmbda*Pival)*Dw))
% cvx_end;

% if strcmp(cvx_status,'Solved')
%     expr=Bw+(-L-L*Pival+Lmbda*Pival)*Dw;
%     Lmbda
%     expr
% end
Cz=1;
cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,p); % n,p
variable Plambda(n,p); % n,p
variable gmma;
minimize(gmma+norm(P*Bw-Y*Dw-Y*Pival*Dw+Plambda*Pival*Dw));
P>=1e-5*eye(n);
[A'*P+P*A-Y*C-(Y*C)',P*Bw-Y*Dw,Cz';
(P*Bw-Y*Dw)',-gmma*eye(1),zeros(1,1);
Cz,zeros(1,1),-gmma*eye(1)]<=0;
cvx_end;

L=inv(P)*Y;
Lmbda=inv(P)*Plambda;
disp(eig(A-L*C))
disp(L)
disp(Lmbda)