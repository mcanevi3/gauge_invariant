%% use PFM to estimate disturbance then feed to Kalman state estimation
% clear; clc; close all;

% 1. Load System Matrices
init; 

% 2. Simulation Setup
N = 1000;              % Number of simulation steps
t = (0:N-1) * tau;     % Time vector
u = sin(2*pi*5*t);     % Example control input

% 3. Noise Covariances
Q_x = diag([1e-6, 1e-6]); 
Q_w = 1e-4;               % Process noise covariance for disturbance
R   = diag([1e-4, 1e-4, 1e-4]); 

% 4. Projection Matrix Setup
Pi = eye(n_y) - C * pinv(C'*C) * C';
Pi_Dw = Pi * Dw;
R_bar = Pi * R * Pi';     % Projected measurement noise covariance

% 5. Full Pre-allocation and Signal Generation
x_true = zeros(n_x, N);
y_meas = zeros(n_y, N);
w_true = zeros(n_w, N);

x_hat  = zeros(n_x, N);
w_hat  = zeros(n_w, N);

% Generate a complex disturbance profile offline
w_true(1, 201:400) = 2.0;                                     % Step
w_true(1, 401:700) = 2.0 + 1.5 * sin(2*pi*2*t(401:700));      % Sine wave offset
w_true(1, 701:N)   = linspace(2.0, -1.0, N-700);              % Ramp down

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
    
    
    % --- STAGE 1: DISTURBANCE ESTIMATION ---
    % Prediction
    w_pred = w_hat(:,k-1);
    P_w_pred = P_w + Q_w;
    
    % Projected measurement
    y_bar = Pi * y_meas(:,k);
    
    % Update
    S_w = Pi_Dw * P_w_pred * Pi_Dw' + R_bar;
    K_w = P_w_pred * Pi_Dw' * pinv(S_w); 
    
    w_hat(:,k) = w_pred + K_w * (y_bar - Pi_Dw * w_pred);
    P_w = (eye(n_w) - K_w * Pi_Dw) * P_w_pred;


    % --- STAGE 2: STATE ESTIMATION ---
    % Prediction (uses previous disturbance estimate w_{k-1|k-1})
    x_pred = Ad * x_hat(:,k-1) + Bud * u(k-1) + Bwd * w_hat(:,k-1);
    P_x_pred = Ad * P_x * Ad' + Q_x;
    
    % Update (uses current disturbance estimate w_{k|k} to correct the innovation)
    S_x = C * P_x_pred * C' + R;
    K_x = P_x_pred * C' / S_x; 
    
    x_hat(:,k) = x_pred + K_x * (y_meas(:,k) - C * x_pred - Dw * w_hat(:,k));
    P_x = (eye(n_x) - K_x * C) * P_x_pred;
end

% 7. Plot Results
figure('Name', 'Two-Stage Kalman Filter Results', 'Position', [100, 100, 800, 800]);

subplot(3,1,1);
plot(t, x_true(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, x_hat(1,:), 'r--', 'LineWidth', 1.5);
ylabel('x_1'); title('State 1: True vs Estimate');
legend('True', 'Estimate', 'Location', 'best');
grid on;

subplot(3,1,2);
plot(t, x_true(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, x_hat(2,:), 'r--', 'LineWidth', 1.5);
ylabel('x_2'); title('State 2: True vs Estimate');
grid on;

subplot(3,1,3);
plot(t, w_true(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(t, w_hat(1,:), 'r--', 'LineWidth', 1.5);
ylabel('w'); xlabel('Time (s)'); title('Disturbance: True vs Estimate');
grid on;