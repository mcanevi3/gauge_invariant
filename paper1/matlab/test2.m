clear;clc;
%%% solving for L then solving for Lmbda
A=-1;
Bw=1;
C=[1;2];
Dw=[1;2]*1+[2;-1]*1;
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

Cz=1;
cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,p); % n,p
variable gmma;
minimize(gmma);
P>=1e-5*eye(n);
[A'*P+P*A-Y*C-(Y*C)',P*Bw-Y*Dw,Cz';
(P*Bw-Y*Dw)',-gmma*eye(1),zeros(1,1);
Cz,zeros(1,1),-gmma*eye(1)]<=0;
cvx_end;

L=inv(P)*Y;
disp("Eig:");
disp(eig(A-L*C))
disp("L:");
disp(L)

%  L=[1,1];


syms x1 x2 x3 x4 x5 x6 real;
Lmbda=[x1,1];
expr=Bw-L*Dw+L*Pival*Dw-Lmbda*Pival*Dw;
x1val=double(solve(expr==0));
disp("Prob:");
disp(vpa(expr,4));
Lmbda=[x1val,1];
expr=Bw-L*Dw+L*Pival*Dw-Lmbda*Pival*Dw;
disp(vpa(expr,4));

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


