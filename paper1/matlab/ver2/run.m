clear;clc;

A=-1;
Bw=1;
C=[1;2];
Dw=[3;1]
% Dw=[1;2]*1+[2;-1]*1;
% C=[1;2;3];
% Bw=[1,3];
% Dw=[[1;2;3]*1+1,[
%    -0.5345   -0.8018
%     0.7745   -0.3382
%    -0.3382    0.4927]*[1;2]];

% % A=[1,0;0,-2];
% % Bw=[1;0];
% % C=[1,0;1,0;0,1e-3];
% % Dw=[1;0;1];

% r=size(Bw,2);
% n=size(A,1);
% p=size(C,1);

% Pival=eye(p)-C*pinv(C'*C)*C';

% angle=acosd(Dw'*C/(norm(C)*norm(Dw)));
% fprintf("Angle: %f\n",angle);

% Cz=eye(n);
% cvx_clear;
% cvx_begin sdp quiet;
% variable P(n,n) symmetric;
% variable Y(n,p); % n,p
% variable gmma;
% minimize(gmma);
% P>=1e-5*eye(n);
% [A'*P+P*A-Y*C-(Y*C)',P*Bw-Y*Dw,Cz';
% (P*Bw-Y*Dw)',-gmma*eye(r),zeros(r,n);
% Cz,zeros(r,n)',-gmma*eye(n)]<=0;
% cvx_end;

% L=inv(P)*Y;
% disp("Eig:");
% disp(eig(A-L*C))
% disp("L:");
% disp(L)
% disp("Bw-LD_w");
% disp(Bw-L*Dw);

% if rank(Pival*Dw)==r
%     disp("rank(Pival*Dw) == r satisfed")
% else
%     disp("rank(Pival*Dw) == r failed")
% end


