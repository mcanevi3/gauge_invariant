clear;clc;
proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
% 1. Define Continuous-Time System Matrices
Ac = [1, 2; -3, -4];
Buc = [1; 0];
Bwc = [0.1; 0.1]*0;
C = [1, 2; 1, 0; 1, 0];
Du = [0; 0; 0];
Dw = [1; 2; 3]; 

% 2. Discretization
Ts = 0.01;
sysC = ss(Ac, [Buc Bwc], C, [Du Dw]); 
sysD = c2d(sysC, Ts); 
Ad = sysD.A; Bud = sysD.B(:, 1); Bwd = sysD.B(:, 2);

% --- NEW: 3. The Projection Matrix ---
% Create the projection matrix onto the left null space of Dw
Pi_Dw = proj(Dw);

% Calculate the virtual measurement matrices
C_star = Pi_Dw * C;

% 4. Tuning Covariances
W = 1.0;       % Tuning weight for state drift caused by unmeasured w
V_noise = 0.05; % Variance of independent sensor noise

Q = Bwd * W * Bwd'; 
% Project the sensor noise covariance into the virtual space
R_star = Pi_Dw * (V_noise * eye(3)) * Pi_Dw'; 

% 5. Initialization
steps = 300; 
x_true = [0; 0];
x_est  = [0; 0];
P      = eye(2);
u      = zeros(1, steps);

x_true_hist = zeros(2, steps);
x_est_hist  = zeros(2, steps);

t = (0:steps-1) * Ts;
w_deterministic = 2 * sin(2 * pi * 1.0 * t); % The road bump

% 6. Projected Kalman Filter Loop
for k = 1:steps
    % --- A. SYSTEM SIMULATION ---
    w_k = w_deterministic(k);          
    v_k = sqrt(V_noise) * randn(3, 1); 
    
    % True measurement (corrupted by massive disturbance)
    y_k = C * x_true + Dw * w_k + v_k;
    
    x_true_hist(:, k) = x_true;
    x_est_hist(:, k)  = x_est; 
    
    % --- B. PROJECT THE MEASUREMENT ---
    % Erase the Dw * w_k term entirely
    y_star = Pi_Dw * y_k; 
    
    % --- C. KALMAN PREDICT ---
    x_predict = Ad * x_est + Bud * u(k);
    P_predict = Ad * P * Ad' + Q;
    
    % --- D. KALMAN CORRECT (Using virtual matrices) ---
    innovation = y_star - C_star * x_predict;
    S_star = C_star * P_predict * C_star' + R_star;
    
    % CRITICAL: Use pinv() because S_star is rank 2 in a 3x3 space
    K = P_predict * C_star' * pinv(S_star); 
    
    x_est = x_predict + K * innovation;
    
    % Update covariance using the Joseph form for numerical stability
    % with reduced-rank matrices
    I_KC = eye(2) - K * C_star;
    P = I_KC * P_predict * I_KC' + K * R_star * K';
    
    % --- E. ADVANCE TRUE SYSTEM ---
    x_true = Ad * x_true + Bud * u(k) + Bwd * w_k;
end

% 7. Plot the Results
figure(1);
subplot(2,1,1);
plot(t, x_true_hist(1,:), 'k', t, x_est_hist(1,:), 'r--', 'LineWidth', 1.5);
legend('True x_1', 'Estimated x_1'); ylabel('State x_1');
title('Kalman Filter with Measurement Disturbance Projection (\Pi_{Dw})'); grid on;

subplot(2,1,2);
plot(t, x_true_hist(2,:), 'k', t, x_est_hist(2,:), 'r--', 'LineWidth', 1.5);
legend('True x_2', 'Estimated x_2'); xlabel('Time (s)'); ylabel('State x_2'); grid on;