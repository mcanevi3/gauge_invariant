clear;clc;

A=[1,2;-3,-4];
B=[1;0];
C=[1,0];
n=size(A,1);note
p = size(C, 1);

Q=eye(n);
Q_sqrt = sqrtm(Q);
R=1;
R_sqrt = sqrtm(R);

Bn = [1;1;1]; 
Bb = [zeros(n, p); -eye(p)]; 
Dn = [zeros(n, p); R_sqrt];

cvx_begin sdp
variable P(n+p, n+p) symmetric
variable Y(n, p)
variable G(p, p)
variable gamma_sq
minimize(gamma_sq)
P >= 1e-6 * eye(n+p); 
AeP = [A*P(1:n, 1:n)-Y*C,  Y;
       G*C,                G];
Caug = [Q_sqrt , zeros(2,1);
        R_sqrt * C , R_sqrt];
[ AeP + AeP',         P*Bn,         P*Bb,         (Caug*P)';
 (P*Bn)',       -gamma_sq*eye(1), zeros(1, p),  Dn'; 
 (P*Bb)',        zeros(p, 1), -gamma_sq*eye(p), zeros(p, p+n); 
  Caug*P,       Dn,           zeros(p+n, p), -gamma_sq*eye(p+n) ] < 0;
cvx_end

if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')
    L = inv(P(1:n, 1:n))*Y;
    Gamma = inv(P(n+1:end, n+1:end))*G;

    disp("L");
    disp(L);
    disp("eig(A-LC):");
    disp(eig(A-L*C));
    disp("Gamma");
    disp(Gamma);
end

