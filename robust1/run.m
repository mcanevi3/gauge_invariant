clear;clc;

A=-1;
Bw=1;
Dw=[1;2]*1+[2;-1]*1;

C0=[1;2];
Ctheta=[1;0];
theta_bar=1;


theta_test=0.5;
C=C0+(theta_test)*Ctheta;

r=size(Bw,2);
n=size(A,1);
p=size(C,1);

Pi0=eye(p)-C0*pinv(C0'*C0)*C0';
Piub=eye(p)-(C0+theta_bar*Ctheta)*pinv((C0+theta_bar*Ctheta)'*(C0+theta_bar*Ctheta))*(C0+theta_bar*Ctheta)';
Pilb=eye(p)-(C0-theta_bar*Ctheta)*pinv((C0-theta_bar*Ctheta)'*(C0-theta_bar*Ctheta))*(C0-theta_bar*Ctheta)';

angle=acosd(Dw'*C/(norm(C)*norm(Dw)));
fprintf("Angle: %f\n",angle);

Cz=eye(n);
cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable S(n,n) symmetric
variable Y(n,p); % n,p
variable gmma;
minimize(gmma);
P>=1e-5*eye(n);
S >= 1e-5*eye(n);
DeltaM = Y*Ctheta;
[ A'*P+P*A - Y*C0 - (Y*C0)' + theta_bar^2*S, ...
  P*Bw - Y*Dw, Cz', DeltaM;
  (P*Bw - Y*Dw)', -gmma*eye(r), zeros(r,n), zeros(r,n);
  Cz, zeros(n,r), -gmma*eye(n), zeros(n,n);
  DeltaM', zeros(n,r), zeros(n,n), -S ] <= 0;
cvx_end;

L=inv(P)*Y;
disp("Eig:");
disp(eig(A-L*C))
disp("L:");
disp(L);

theta_vals = linspace(-theta_bar,theta_bar,10);
eig_vals = zeros(1,length(theta_vals));

for k = 1:length(theta_vals)
    Ck = C0 + theta_vals(k)*Ctheta;
    eig_vals(k) = eig(A - L*Ck);   % L is fixed
end
disp("Eig_vals:");
disp(eig_vals);

% disp("Bw-LD_w");
% disp(Bw-L*Dw);

% if rank(Pival*Dw)==r
%     disp("rank(Pival*Dw) == r satisfed")
% else
%     disp("rank(Pival*Dw) == r failed")
% end


