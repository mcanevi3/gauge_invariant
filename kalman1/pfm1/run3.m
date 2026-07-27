%% PFM kalman for disturbance, PFM kalman for state
clear; clc; close all;
rng(42, 'twister');
% 1. Load System Matrices
init; 

% --- For Disturbance Estimator (Nullifies C) ---
Pi_C = eye(n_y) - C * pinv(C'*C) * C';
Pi_C_Dw = Pi_C * Dw;
R_bar_w = Pi_C * R * Pi_C';

% --- For State Estimator (Nullifies Dw) ---
Pi_D = eye(n_y) - Dw * pinv(Dw'*Dw) * Dw';
Pi_D_C = Pi_D * C;
R_bar_x = Pi_D * R * Pi_D';

% 5. Full Pre-allocation and Signal Generation
x_true = zeros(n_x, N);
y_meas = zeros(n_y, N);


x_hat  = zeros(n_x, N);
w_hat  = zeros(n_w, N);

% % Generate a complex disturbance profile offline
% w_true(1, 201:400) = 2.0;                                     
% w_true(1, 401:700) = 2.0 + 1.5 * sin(2*pi*2*t(401:700));      
% w_true(1, 701:N)   = linspace(2.0, -1.0, N-700);              

% Pre-generate noise arrays using element-wise multiplication (.*)
V_x = sqrt(diag(Q_x)) .* randn(n_x, N);
V_k = sqrt(diag(R)) .* randn(n_y, N);

% Initialize states and covariances
x_true(:,1) = [0; 0];
x_hat(:,1)  = [0; 0];
w_hat(:,1)  = 0;

P_x = eye(n_x);
P_w = eye(n_w);

% 6. Main Simulation and Filter Loop
for k = 2:N
    % --- SIMULATE TRUE SYSTEM ---
    x_true(:,k) = Ad * x_true(:,k-1) + Bud * u(k-1) + Bwd * w_true(:,k-1) + V_x(:,k);
    y_meas(:,k) = C * x_true(:,k) + Dw * w_true(:,k) + V_k(:,k);
    
    
    % --- ESTIMATOR 1: DISTURBANCE (PFM nullifying C) ---
    w_pred = w_hat(:,k-1);
    P_w_pred = P_w + Q_w;
    
    y_bar_w = Pi_C * y_meas(:,k);
    
    S_w = Pi_C_Dw * P_w_pred * Pi_C_Dw' + R_bar_w;
    K_w = P_w_pred * Pi_C_Dw' * pinv(S_w); 
    
    w_hat(:,k) = w_pred + K_w * (y_bar_w - Pi_C_Dw * w_pred);
    P_w = (eye(n_w) - K_w * Pi_C_Dw) * P_w_pred;


    % --- ESTIMATOR 2: STATE (PFM nullifying Dw) ---
    % Prediction (still requires w_hat for the physical system propagation)
    x_pred = Ad * x_hat(:,k-1) + Bud * u(k-1) + Bwd * w_hat(:,k-1);
    P_x_pred = Ad * P_x * Ad' + Q_x;
    
    % Projected Measurement for State
    y_bar_x = Pi_D * y_meas(:,k);
    
    % Update (Innovation is now structurally immune to Dw * w_k)
    S_x = Pi_D_C * P_x_pred * Pi_D_C' + R_bar_x;
    K_x = P_x_pred * Pi_D_C' * pinv(S_x); % pinv needed due to projection
    
    x_hat(:,k) = x_pred + K_x * (y_bar_x - Pi_D_C * x_pred);
    P_x = (eye(n_x) - K_x * Pi_D_C) * P_x_pred;
end


%% save results
results_uio.t = t;
results_uio.x_true = x_true;
results_uio.w_true = w_true;
results_uio.x_hat = x_hat;
results_uio.w_hat = w_hat;

% Save the struct to a binary .mat file (fast and compressed)
save('sim_results_uio.mat', '-struct', 'results_uio');