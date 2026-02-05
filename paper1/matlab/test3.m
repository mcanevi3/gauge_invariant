clear;clc;
rng(233)
% is there a solution Bw-L*Dw\neq 0 ???

n=2; % num of states
r=1; % size of w
p=3; % num of sensors
A=rand(n,n);
Bw=rand(n,r);
C=rand(p,n);

Cperp=null(C');
Pival=eye(p)-C*pinv(C'*C)*C';

Dw=rand(p,r);
% Dw=C/norm(C)*[1;1]+Cperp/norm(Cperp)*1;
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
disp(L);

expr=Bw-L*Dw;
disp("Bw-L*Dw:");
disp(expr);
