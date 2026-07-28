%% Augmented Kalman estimating w and x at once
clear; clc; close all;

% 1. Lock Random Seed for Reproducibility
rng(43, 'twister'); 

% 2. Load System Matrices
init; 

% --- Augmented System Matrices ---
A_aug = [Ad, Bwd; 
         zeros(n_w, n_x), eye(n_w)];
B_aug = [Bud; 
         zeros(n_w, size(Bud, 2))];
C_aug = [C, Dw];
Q_aug = blkdiag(Q_x, Q_w);

n_z = size(A_aug, 1); % Total augmented states (n_x + n_w)

% Signal Pre-allocation
x_true = zeros(n_x, N);
y_meas = zeros(n_y, N);
z_hat  = zeros(n_z, N); % Stores both x_hat and w_hat


% Pre-generate identical noise arrays for fair comparison
V_x = sqrt(diag(Q_x)) .* randn(n_x, N);
V_k = sqrt(diag(R)) .* randn(n_y, N);

% Initialization
x_true(:,1) = [0; 0];
z_hat(:,1)  = [0; 0; 0]; % Initial augmented state: [x1; x2; w]
P_aug = eye(n_z);

eps_k_hist=zeros(1,N);

% 4. Main Simulation and Filter Loop
for k = 2:N
    % --- SIMULATE TRUE SYSTEM ---
    x_true(:,k) = Ad * x_true(:,k-1) + Bud * u(k-1) + Bwd * w_true(:,k-1) + V_x(:,k);
    y_meas(:,k) = C * x_true(:,k) + Dw * w_true(:,k) + V_k(:,k);
    
    % --- AUGMENTED KALMAN FILTER ---
    % Prediction
    z_pred = A_aug * z_hat(:,k-1) + B_aug * u(k-1);
    P_pred = A_aug * P_aug * A_aug' + Q_aug;
    
    % Update
    S_k = C_aug * P_pred * C_aug' + R;
    K_k = P_pred * C_aug' / S_k; % Standard division is fine here (no projection)
    
    % Normalized Innovation Squared (NIS)
    v_k=y_meas(:,k) - C_aug * z_pred;
    eps_k = v_k'*inv(S_k)*v_k;
    eps_k_hist(k)=eps_k;

    z_hat(:,k) = z_pred + K_k * (y_meas(:,k) - C_aug * z_pred);
    P_aug = (eye(n_z) - K_k * C_aug) * P_pred;
end

% 5. Extract and Package Results Efficiently
% Split z_hat back into x_hat and w_hat for easy comparison with the UIO method
x_hat = z_hat(1:n_x, :);
w_hat = z_hat(n_x+1:end, :);

results_aug.t = t;
results_aug.x_true = x_true;
results_aug.w_true = w_true;
results_aug.x_hat = x_hat;
results_aug.w_hat = w_hat;
results_aug.eps_k_hist = eps_k_hist;

% Save the struct to a binary .mat file
save('sim_results_aug.mat', '-struct', 'results_aug');
disp('Simulation complete. Results saved to sim_results_aug.mat');