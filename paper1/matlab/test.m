clear;clc;
%%% L fixed using ackerman minimize norm over Lmbda
A=-1;
Bw=1;
C=[1;2];
Dw=[1;2];
L=[1,1];
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

% syms x1 x2 x3 x4 x5 x6 real;
% Lmbda=[x1,x2,x3;x4,x5,x6];
% expr=Bw+(-L-L*Pival+Lmbda*Pival)*Dw;
% disp("Prob:");
% disp(expr);


cvx_clear;
cvx_begin sdp quiet;
variable Lmbda(n,p)
minimize(norm(Bw+(-L-L*Pival+Lmbda*Pival)*Dw))
cvx_end;

if strcmp(cvx_status,'Solved')
    expr=Bw+(-L-L*Pival+Lmbda*Pival)*Dw;
    Lmbda
    expr
end