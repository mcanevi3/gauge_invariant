%% Kalman filter augmented with disturbance x1,x2,w
clear;clc;
init;

% Augment the state vector z = [x1; x2; w]
A_aug = [Ad, Bwd; 
         zeros(n_w, n_x), eye(n_w)]; % Random walk for w in discrete-time
B_aug = [Bud; 
         zeros(n_w, size(Bud, 2))];
C_aug = [C, Dw];

%% luenberger
% poles_luen = [0.80, 0.81, 0.82]; 
% L_luen=place(A_aug',C_aug',poles_luen)';
% L_luen =[
%    -0.0009    0.0013   -0.0005
%     0.0105   -0.0021   -0.0021
%     0.0701    0.1421    0.2140
% ];
des_luen;

%% kalman 
Q_w_var = 1.0;   % Tuning: Expected rate of change of disturbance
R_v_var = 0.05;  % Tuning: Sensor noise variance
% Process noise applied primarily to the random walk of the disturbance
Q_kf = blkdiag(zeros(n_x), Q_w_var * eye(n_w)); 
R_kf = R_v_var * eye(n_y);

%% Simulation
N_steps = 3000; 
t = (0:N_steps-1) * tau;

u = zeros(1, N_steps);
w_true = 2 * sin(2 * pi * 10.0 * t);% + t.^2; 

x_true_hist = zeros(n_x + n_w, N_steps);
x_luen_hist = zeros(n_x + n_w, N_steps);
x_kf_hist   = zeros(n_x + n_w, N_steps);
% x_pi_hist   = zeros(n_x + n_w, N_steps);

% Initial Conditions
x_true = zeros(n_x, 1);          % True physical states
x_luen = zeros(n_x + n_w, 1);    % Luenberger augmented estimate
x_kf   = zeros(n_x + n_w, 1);    % Kalman augmented estimate
P_kf   = eye(n_x + n_w);         % Initial Kalman covariance

for k = 1:N_steps
    % ==========================================
    % A. TRUE SYSTEM DYNAMICS & MEASUREMENT
    % ==========================================
    % Construct current augmented true state for logging and output
    x_aug_true = [x_true; w_true(k)];
    
    % Generate noisy measurement
    v_k = sqrt(R_v_var) * randn(n_y, 1); 
    y_k = C_aug * x_aug_true + v_k;
    
    % Log current states
    x_true_hist(:, k) = x_aug_true;
    x_luen_hist(:, k) = x_luen;
    x_kf_hist(:, k)   = x_kf;
    
    % ==========================================
    % B. LUENBERGER OBSERVER UPDATE
    % ==========================================
    y_hat_luen = C_aug * x_luen;
    x_luen = A_aug * x_luen + B_aug * u(k) + L_luen * (y_k - y_hat_luen);
    
    % ==========================================
    % C. KALMAN FILTER UPDATE (Predict -> Correct)
    % ==========================================
    % 1. Predict
    x_predict = A_aug * x_kf + B_aug * u(k);
    P_predict = A_aug * P_kf * A_aug' + Q_kf;
    
    % 2. Correct
    S_cov = C_aug * P_predict * C_aug' + R_kf;
    K_gain = P_predict * C_aug' / S_cov; % More robust than inv(S)
    
    x_kf = x_predict + K_gain * (y_k - C_aug * x_predict);
    P_kf = (eye(n_x + n_w) - K_gain * C_aug) * P_predict;
    
    % ==========================================
    % D. ADVANCE TRUE SYSTEM (Physical states only)
    % ==========================================
    x_true = Ad * x_true + Bud * u(k) + Bwd * w_true(k);
end

plot_res;