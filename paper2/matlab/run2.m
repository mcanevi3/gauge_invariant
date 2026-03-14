clear;clc;

A=[1,2;-3,-4];
Bw=[2;1];
C=[1,0;1,0;1,1];
Dw=[1;2;3];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);
Cperp= null(C');

Dw=Cperp*1+0*C*[1;1];

n = size(A,1);
p = size(C,1);

alpha = 2;   % desired stability margin

cvx_begin sdp quiet;

variable P(n,n) symmetric
variables Y(n,p) Z(n,n) gmma

minimize(gmma)

subject to

P >= 1e-6*eye(n);

% Stability region: Re(lambda)<-alpha
A'*P + P*A - C'*Y' - Y*C + 2*alpha*P <= -1e-6*eye(n);

% Disturbance decoupling
Y*Dw == P*Bw;

% Definition of Z
Z == Y*C;

% Norm bound ||Z||_2 <= gamma
[gmma*eye(n) Z;
 Z' gmma*eye(n)] >= 0;

cvx_end

disp("Status:"+string(cvx_status));
% Recover observer gain
L = P\Y;

disp("L:");
disp(L);
disp("LC:");
disp(L*C);
disp("||LC||:"+string(norm(L*C,inf)));
disp("Eig:");
disp(eig(A-L*C));
disp("Bw-L*Dw");
disp(Bw-L*Dw);


