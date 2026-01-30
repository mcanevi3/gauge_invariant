clear;clc;

A=[1,2;-3,-4];
Bw=[1;1];
C=[2,1;4,2];
Dw=[1;1];

Cz=[1,0;0,1];
% Dimensions
n = size(A,1);
p = size(C,1);
r = size(Bw,2);
qz = size(Cz,1);

cvx_begin sdp quiet
    cvx_precision high

    % Decision variables
    variable P1(n,n) symmetric
    variable P2(p,p) symmetric
    variable P12(n,p)
    variable Y(n,p)
    variable Yi(n,p)
    variable gmma

    % Lyapunov matrix
    P = [P1  P12;
         P12' P2];

    % Objective
    minimize(gmma)
    gmma >= 1e-9;

    % Positivity constraints
    P1 >= 1e-9*eye(n);
    P2 >= 1e-9*eye(p);
    P  >= 1e-9*eye(n+p);

    % LMI blocks
    Psi11 = P1*A + A'*P1 ...
            - Y*C - C'*Y' ...
            + P12*C + C'*P12';

    Psi12 = A'*P12 - Yi;

    Psi22 = P12'*C + C'*P12;

    Psi13 = P1*Bw - Y*Dw;
    Psi23 = P12'*Bw + P2*Dw;

    % Bounded real LMI
    LMI = [ Psi11   Psi12    Psi13      Cz';
            Psi12'  Psi22    Psi23      zeros(p,qz);
            Psi13'  Psi23'   -gmma*eye(r)  zeros(r,qz);
            Cz      zeros(qz,p) zeros(qz,r) -gmma*eye(qz) ];

    LMI <= -1e-12*eye(size(LMI));

cvx_end

if strcmp(cvx_status,"Solved") || strcmp(cvx_status,"Inaccurate/Solved")
    L  = P1 \ Y;
    Li = P1 \ Yi;

    Acl=[A-L*C,-Li;C,zeros(2,2)];

else
    disp("Status:"+string(cvx_status))

end