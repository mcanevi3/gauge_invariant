clear;clc;
% 1. Define Continuous-Time System Matrices
Ac = [1, 2; -3, -4];
Buc = [1; 0];
Bwc = [0.1; 0.1];
C = [1, 2; 1, 0; 1, 0];
Du = [0; 0; 0];
Dw = [1; 2; 3]; 

% 2. Discretization using zero-order hold
Ts = 0.001;
sysC = ss(Ac, [Buc Bwc], C, [Du Dw]); 
sysD = c2d(sysC, Ts); 

Ad = sysD.A;
Bud = sysD.B(:, 1);
Bwd = sysD.B(:, 2);

% --- NEW: 3. State Augmentation ---
% Augment the state vector z = [x1; x2; w]
A_aug = [Ad, Bwd; 
         0, 0, 1]; 
B_aug = [Bud; 
         0];
C_aug = [C, Dw];

% 4. Tuning the Augmented Filter
% W_xi tunes how fast we think the unknown disturbance w can change.
% V_noise is the independent sensor static.
W_xi = 1.0;     
V_noise = 0.05; 

% Q_aug now only applies process noise to the 3rd state (the disturbance w)
Q_aug = diag([0, 0, W_xi]); 
R_aug = V_noise * eye(3); 

% 5. Initialization
steps = 3000; 
z_true = [0; 0; 0]; % True augmented state
z_est  = [0; 0; 0]; % Estimated augmented state
P      = eye(3);
u      = zeros(1, steps);

z_true_hist = zeros(3, steps);
z_est_hist  = zeros(3, steps);

t = (0:steps-1) * Ts;
% The true disturbance is a sine wave (but the filter doesn't know this!)
w_deterministic = 2 * sin(2 * pi * 10.0 * t)+t.^2; 

% 6. Augmented Kalman Filter Loop (Predict-Correct Form)
for k = 1:steps
    % --- A. SYSTEM SIMULATION ---
    % Inject the deterministic w into the true state vector for logging
    z_true(3) = w_deterministic(k); 
    v_k = sqrt(V_noise) * randn(3, 1); 
    
    % True measurement (calculated using augmented matrices)
    y_k = C_aug * z_true + v_k;
    
    z_true_hist(:, k) = z_true;
    z_est_hist(:, k)  = z_est; 
    
    % --- B. KALMAN PREDICT ---
    z_predict = A_aug * z_est + B_aug * u(k);
    P_predict = A_aug * P * A_aug' + Q_aug;
    
    % --- C. KALMAN CORRECT ---
    innovation = y_k - C_aug * z_predict;
    S = C_aug * P_predict * C_aug' + R_aug;
    K = P_predict * C_aug' * inv(S); 
    
    z_est = z_predict + K * innovation;
    P     = (eye(3) - K * C_aug) * P_predict;
    
    % --- D. ADVANCE TRUE SYSTEM (Only the physical states x1, x2) ---
    z_true(1:2) = Ad * z_true(1:2) + Bud * u(k) + Bwd * w_deterministic(k);
end

% 7. Plot the Results
figure(1);
subplot(3,1,1);cla;hold on;grid on;
xlabel('t');ylabel('x_1');title('Augmented Filter: State and Disturbance Estimation'); 
plot(t, z_true_hist(1,:), 'k', 'LineWidth', 3);
plot(t, z_est_hist(1,:), 'r', 'LineWidth', 2);
legend('True x_1', 'Estimated x_1'); 

subplot(3,1,2);cla;hold on;
xlabel('t');ylabel('x_2'); grid on;
plot(t, z_true_hist(2,:), 'k', 'LineWidth', 3);
plot(t, z_est_hist(2,:), 'r', 'LineWidth', 2);
legend('True x_2', 'Estimated x_2'); 

subplot(3,1,3);cla;hold on;grid on;
xlabel('t'); ylabel('w');
plot(t, z_true_hist(3,:), 'k', 'LineWidth', 3);
plot(t, z_est_hist(3,:), 'b', 'LineWidth', 2);
legend('True Disturbance w', 'Estimated Disturbance w');
