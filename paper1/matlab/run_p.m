clear;clc;

A=[1,2;-3,-4];
Bw=[1;1];
C=[2,1;4,2];
Dw=[1;1];

Cz=[1,0;0,1];

% Dimensions
n  = size(A,1);      % state dimension
p  = size(C,1);      % output dimension
r  = size(Bw,2);     % disturbance dimension
qz = size(Cz,1);     % performance output dimension

cvx_begin sdp quiet
    cvx_precision high

    % Decision variables
    variable P(n,n) symmetric
    variable Y(n,p)
    variable gmma

    % Objective
    minimize(gmma)

    % Positivity
    P >= 1e-6*eye(n);
    gmma >= 1e-6;

    % LMI blocks
    Phi = P*A + A'*P - Y*C - C'*Y';

    % Bounded real LMI
    LMI = [ Phi,            P*Bw - Y*Dw,   Cz';
            (P*Bw - Y*Dw)', -gmma*eye(r), zeros(r,qz);
            Cz,             zeros(qz,r),   -gmma*eye(qz) ];

    LMI <= -1e-6*eye(size(LMI));

cvx_end


if strcmp(cvx_status,"Solved") || strcmp(cvx_status,"Inaccurate/Solved")
    L = P \ Y;
else
    disp("Status:"+string(cvx_status))
end
