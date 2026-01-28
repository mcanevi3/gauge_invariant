clear; clc; close all;

%% ---------------- USER: Define plant ----------------
A = [0 1; -1 -0.3];
B = [0; 1];
C = [1 0];

n = size(A,1);
m = size(B,2);
p = size(C,1);

%% ---------------- Bias model ----------------
bias_mode = "drift";   % "constant" or "drift"
b0 = 0.20*ones(p,1);

%% ---------------- Simulation horizon ----------------
Tend = 10;
dt_noise = 1e-3;
tgrid = (0:dt_noise:Tend)';

%% ---------------- Disturbance + noise base levels ----------------
sigma_w_base  = 1e-2;    % process noise std
sigma_n_base  = 5e-2;    % measurement noise std
sigma_wb_base = 5e-4;    % bias drift std (if drift mode)

%% ---------------- Input (control) ----------------
u_fun = @(t, x) zeros(m,1); % default: u=0

%% ---------------- Classical observer gain (baseline) ----------------
% Use pole placement just for baseline comparison.
poles_classical = -[2 3];
Lc = place(A', C', poles_classical)';

%% ---------------- Augmented model ----------------
Aa = [A zeros(n,p);
      zeros(p,n) zeros(p,p)];
Ba = [B;
      zeros(p,m)];
Ca = [C eye(p)];
na = n + p; % augmented dimension

%% ---------------- Choose performance output z = Ce * e_a ----------------
% You can weight bias error less/more by changing wbias.
wbias = 0.5; % bias error weight (try 0.1, 0.5, 1.0)
Ce = blkdiag(eye(n), wbias*eye(p));   % z = [e_x; wbias*e_b]

% Disturbance input to error system: e_a_dot = (Aa - La Ca)e_a + Be w
% Use conservative Be = I (generic additive disturbance).
Be = eye(na);

%% ---------------- CVX DESIGN: minimize gamma ----------------
fprintf("\n--- CVX: Designing GI observer gain La via SDP/LMI ---\n");
[La, gamma_opt, Popt] = gi_observer_cvx_design(Aa, Ca, Be, Ce);

fprintf("CVX done. Optimal gamma ~ %.6g\n", gamma_opt);
disp("GI observer gain La (augmented) = "); disp(La);

%% ---------------- Stress-test sweep ----------------
alphas = 0:0.25:1; %5
rms_ex_class = zeros(size(alphas));
rms_ex_gi    = zeros(size(alphas));
rms_eb_gi    = zeros(size(alphas));

rng(7);

for k = 1:numel(alphas)
    alpha = alphas(k);

    W  = sigma_w_base*alpha  * randn(numel(tgrid), n);
    N  = sigma_n_base*alpha  * randn(numel(tgrid), p);
    Wb = sigma_wb_base*alpha * randn(numel(tgrid), p);

    x0 = zeros(n,1);
    xhat0 = zeros(n,1);
    xahat0 = [zeros(n,1); zeros(p,1)];
    b_init = b0;

    P = struct();
    P.A = A; P.B = B; P.C = C;
    P.Aa = Aa; P.Ba = Ba; P.Ca = Ca;
    P.Lc = Lc; P.La = La;
    P.u_fun = u_fun;
    P.bias_mode = bias_mode;

    P.tgrid = tgrid;
    P.W = W; P.N = N; P.Wb = Wb;

    s0 = [x0;
          xhat0;
          xahat0;
          b_init];

    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    [t, s] = ode45(@(t,s) odefun_all(t,s,P), [0 Tend], s0, opts);

    x      = s(:, 1:n);
    xhat_c = s(:, n+1:2*n);
    xahat  = s(:, 2*n+1:2*n+(n+p));
    btraj  = s(:, 2*n+(n+p)+1:end);

    xhat_gi = xahat(:, 1:n);
    bhat_gi = xahat(:, n+1:end);

    ex_c = x - xhat_c;
    ex_g = x - xhat_gi;
    eb_g = btraj - bhat_gi;

    rms_ex_class(k) = sqrt(mean(sum(ex_c.^2,2)));
    rms_ex_gi(k)    = sqrt(mean(sum(ex_g.^2,2)));
    rms_eb_gi(k)    = sqrt(mean(sum(eb_g.^2,2)));

    fprintf("alpha=%.2f | RMS(x-xhat): class=%.4g, GI=%.4g | RMS(b-bhat)=%.4g\n", ...
        alpha, rms_ex_class(k), rms_ex_gi(k), rms_eb_gi(k));
end

%% ---------------- Plots ----------------
figure; plot(alphas, rms_ex_class, '-o', alphas, rms_ex_gi, '-s', 'LineWidth', 1.2);
grid on; xlabel('\alpha (noise scale)'); ylabel('RMS state estimation error');
legend('Classical observer','GI observer (CVX-LMI)','Location','northwest');
title('Noise robustness: RMS state estimation error');

figure; plot(alphas, rms_eb_gi, '-d', 'LineWidth', 1.2);
grid on; xlabel('\alpha (noise scale)'); ylabel('RMS bias estimation error');
title('GI observer (CVX-LMI): bias estimation error vs noise scale');

%% ======================= Local functions =======================

function [La, gamma_opt, Popt] = gi_observer_cvx_design(Aa, Ca, Be, Ce)
    % Designs La via CVX by minimizing gamma subject to BRL-LMI:
    % [Acl'P + P Acl,  PBe,  Ce';
    %  Be'P,          -gI,   0;
    %  Ce,            0,    -gI] < 0
    %
    % Convexify with Y = P La:
    % Acl'P + P Acl = Aa'P + P Aa - Ca'Y' - Y Ca
    %
    % Output:
    %   La = P\Y
    %   gamma_opt = optimal gamma
    %   Popt = optimal P

    na = size(Aa,1);
    p  = size(Ca,1);
    nw = size(Be,2);
    nz = size(Ce,1);

    epsPD = 1e-7;      % small positive number for strictness
    epsLMI = 1e-7;

    cvx_begin sdp quiet
        variable P(na,na) symmetric
        variable Y(na,p)
        variable gm(1,1)

        minimize( gm )
        
        subject to
%            gm >= epsPD;
            gm>=2;
            P >= epsPD*eye(na);

            % LMI block
            Xi = Aa'*P + P*Aa - Ca'*Y' - Y*Ca;

            M = [ Xi,        P*Be,         Ce';
                  Be'*P,     -gm*eye(nw), zeros(nw,nz);
                  Ce,        zeros(nz,nw),  -gm*eye(nz) ];

            M <= -epsLMI*eye(size(M,1));
    cvx_end

    if ~strcmpi(cvx_status, 'Solved') && ~strcmpi(cvx_status, 'Inaccurate/Solved')
        error("CVX did not solve the SDP. Status: %s", cvx_status);
    end

    Popt = P;
    gamma_opt = gm;
    La = P \ Y; % La = inv(P)*Y (more numerically stable)
end

function ds = odefun_all(t, s, P)
    n = size(P.A,1);
    p = size(P.C,1);

    idx = 0;
    x = s(idx+1:idx+n); idx = idx+n;
    xhat_c = s(idx+1:idx+n); idx = idx+n;
    xahat = s(idx+1:idx+(n+p)); idx = idx+(n+p);
    b = s(idx+1:idx+p);

    w  = interp_noise(t, P.tgrid, P.W);   % n x 1
    nmeas = interp_noise(t, P.tgrid, P.N);% p x 1
    wb = interp_noise(t, P.tgrid, P.Wb);  % p x 1

    u = P.u_fun(t, x);

    xdot = P.A*x + P.B*u + w;

    if P.bias_mode == "constant"
        bdot = zeros(p,1);
    else
        bdot = wb;
    end

    y = P.C*x + b + nmeas;

    xhat_c_dot = P.A*xhat_c + P.B*u + P.Lc*(y - P.C*xhat_c);

    xahat_dot = P.Aa*xahat + P.Ba*u + P.La*(y - P.Ca*xahat);

    ds = [xdot;
          xhat_c_dot;
          xahat_dot;
          bdot];
end

function y = interp_noise(t, tgrid, M)
    y = interp1(tgrid, M, t, 'linear', 'extrap')';
end

