clear;clc;

A=[1,2;-3,-4];
B=[1;0];
C=[1,0];
n=size(A,1);

Q=eye(2);
R=1;

cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable Y(n,1);
variable gm;
minimize(gm);
gm>=0;
P>=0;
[(A*P-Y*C)+(A*P-Y*C)'+Q,P,(C*P)'*sqrtm(R);
  P,-gm*eye(n),zeros(n,1);
  sqrtm(R)*(C*P),zeros(1,n),-gm*eye(1)]<=0;
cvx_end;
if strcmp(cvx_status,'Solved')
    L=inv(P)*Y;
    
    disp("L");
    disp(L);
    
    disp("eig(A-LC):");
    disp(eig(A-L*C));

end

