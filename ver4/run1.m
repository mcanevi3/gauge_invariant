clear;clc;

% A=-1;
% Bw=1;
% C=[2;1];

% A=[1,2;-3,-4];
% Bw=[1,1;1,1];
% C=[2,1;1,2];

%-% A=[1,2;-3,-4];
%-% Bw=eye(2);
%-% C=[2;1];
proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';



A=[0,0,1,0;0,0,0,1;31.4752,3.0426,-0.4107,96.7383;-13.6919,9.4947,0.3349,165.0421];
C=[1,0,0,0;0,1,0,0;0,0,1,0];
Bw=[1;1;1];
Dw=[1;1;1];

n=size(A,1);
p=size(C,1);

cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,p); % n,p
variable gmma;
minimize(gmma);
P>=1e-3*eye(n);
[A'*P+P*A-Y*C-C'*Y',P*Bw-Y*Dw,Cz';
(P*Bw-Y*Dw)',-gmma*eye(p),zeros(p,p);
Cz,zeros(p,p),-gmma*eye(p)]<=0;
cvx_end;

L=inv(P)*Y;
gmma=sqrt(ro);
disp(eig(A-L*C))
G=ss(A-L*C,(Bw-L)*C,eye(n),zeros(n,n));
L

Pi=eye(2)-C'*C/(C*C')