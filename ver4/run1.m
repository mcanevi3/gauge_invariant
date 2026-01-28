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

n=size(A,1);
p=size(C,1);

cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,p); % n,p
variable ro;
minimize(ro);
P>=1e-3*eye(n);
[A'*P+P*A-Y*C-C'*Y',(P*Bw-Y)*C;C'*(P*Bw-Y)',-ro*eye(2)]<=0;
cvx_end;

L=inv(P)*Y;
gmma=sqrt(ro);
disp(eig(A-L*C))
G=ss(A-L*C,(Bw-L)*C,eye(n),zeros(n,n));
L

Pi=eye(2)-C'*C/(C*C')