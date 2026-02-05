clear;clc;
%%% solving for L then solving for Lmbda
A=-1;
Bw=1;
C=[1;2];
Dw=[1;2]*1+[2;-1]*1;
Dw=[1;2]*1+[2;-1]*0.2;
% we need a orthogonal component if its aligned with C it perishes
% Dw=[0.4472;0.8944];

% A=[1,0;0,-2];
% Bw=[1;0];
% C=[1,0;1,0;0,1e-3];
% Dw=[1;0;1];


r=size(Bw,2);
n=size(A,1);
p=size(C,1);

Pival=eye(p)-C*pinv(C'*C)*C';
% Dw=null(Pival);

angle=acosd(Dw'*C/(norm(C)*norm(Dw)));
fprintf("Angle: %f\n",angle);

Cz=eye(n);
cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,p); % n,p
variable gmma;
minimize(gmma);
P>=1e-5*eye(n);
[A'*P+P*A-Y*C-(Y*C)',P*Bw-Y*Dw,Cz';
(P*Bw-Y*Dw)',-gmma*eye(r),zeros(r,n);
Cz,zeros(r,n)',-gmma*eye(n)]<=0;
cvx_end;

L=inv(P)*Y;
disp("Eig:");
disp(eig(A-L*C))
disp("L:");
disp(L)
disp("Bw-LD_w");
disp(Bw-L*Dw)
%  L=[1,1];

rank(Pival*Dw)

% syms x1 x2 x3 x4 x5 x6 real;
% Lmbda=[x1,1];
% expr=Bw-L*Dw+L*Pival*Dw-Lmbda*Pival*Dw;
% x1val=double(solve(expr==0));
% disp("Prob:");
% disp(vpa(expr,4));
% Lmbda=[x1val,1];
% expr=Bw-L*Dw+L*Pival*Dw-Lmbda*Pival*Dw;
% disp(vpa(expr,4));

% Lmbda=Lmbda*0;
% Lmbda=[x1val,1]*0;
% expr=Bw+(-L-L*Pival+Lmbda*Pival)*Dw

% cvx_clear;
% cvx_begin sdp quiet;
% variable Lmbda(n,p)
% minimize(norm(Bw-L*Dw+(L-Lmbda)*Dw))
% cvx_end;

% if strcmp(cvx_status,'Solved')
%     expr=Bw-L*Dw+(L-Lmbda)*Dw;
%     expr
% end


