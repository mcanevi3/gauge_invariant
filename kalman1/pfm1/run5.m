%% Adaptive Augmented Kalman estimating w and x at once adapting Sk
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

Q_x = diag([1e-6, 1e-6]); 
Q_w = 1e-4;               
R   = diag([1e-4, 1e-4, 1e-4]); 
Q_aug = blkdiag(Q_x, Q_w);

n_z = size(A_aug, 1); % Total augmented states (n_x + n_w)

% % 3. Simulation Parameters
% N = 1000;              
% t = (0:N-1) * tau;     
% u = sin(2*pi*5*t); 

% Signal Pre-allocation
x_true = zeros(n_x, N);
y_meas = zeros(n_y, N);
% w_true = zeros(n_w, N);
z_hat  = zeros(n_z, N); % Stores both x_hat and w_hat

% % Malicious Disturbance Profile (to trigger the adaptation)
% w_true(1, 201:500) = 5 * sign(sin(2*pi*50*t(201:500))); % High-freq chatter
% w_true(1, 600) = 15;  % Massive shock
% w_true(1, 700) = -20;
% w_true(1, 800) = 25;
% w_true(1, 900) = -15;

% Pre-generate identical noise arrays
V_x = sqrt(diag(Q_x)) .* randn(n_x, N);
V_k = sqrt(diag(R)) .* randn(n_y, N);

% Initialization
x_true(:,1) = [0; 0];
z_hat(:,1)  = [0; 0; 0]; 
P_aug = eye(n_z);
eps_k_hist = zeros(1, N);

% --- Adaptive Tuning Parameters ---
% For n_y = 3 measurements, the 95% confidence interval for a Chi-squared 
% distribution is approx 7.81. This is our fault detection threshold.
tau_chi2 = 7.81; 
inflation_scale = 100; % Tuning gain for how aggressively R inflates

% 4. Main Simulation and Filter Loop
for k = 2:N
    % --- SIMULATE TRUE SYSTEM ---
    x_true(:,k) = Ad * x_true(:,k-1) + Bud * u(k-1) + Bwd * w_true(:,k-1) + V_x(:,k);
    y_meas(:,k) = C * x_true(:,k) + Dw * w_true(:,k) + V_k(:,k);
    
    % --- ADAPTIVE AUGMENTED KALMAN FILTER ---
    % 1. Prediction
    z_pred = A_aug * z_hat(:,k-1) + B_aug * u(k-1);
    P_pred = A_aug * P_aug * A_aug' + Q_aug;
    
    % 2. PASS 1: Nominal NIS Calculation
    S_nom = C_aug * P_pred * C_aug' + R;
    v_k = y_meas(:,k) - C_aug * z_pred;
    eps_k = v_k' * inv(S_nom) * v_k;
    eps_k_hist(k) = eps_k;
    
    % 3. PASS 2: Adaptive Covariance Inflation
    if eps_k > tau_chi2
        % Fault detected! Inflate R strictly along the Dw vector geometry.
        % The further NIS is above the threshold, the harder we project it out.
        gamma_k = inflation_scale * (eps_k - tau_chi2);
        R_adapt = R + gamma_k * (Dw * Dw');
        
        % Recompute innovation covariance with the inflated R
        S_k = C_aug * P_pred * C_aug' + R_adapt;
    else
        % Normal Gaussian operation
        S_k = S_nom;
    end

    % 4. Final Update
    K_k = P_pred * C_aug' / S_k; 
    z_hat(:,k) = z_pred + K_k * v_k;
    
    % Note: Using the standard Joseph form or simple form is acceptable here, 
    % but ensure it uses the final K_k.
    P_aug = (eye(n_z) - K_k * C_aug) * P_pred;
end

% 5. Extract and Package Results 
x_hat = z_hat(1:n_x, :);
w_hat = z_hat(n_x+1:end, :);

results_aug_adapt.t = t;
results_aug_adapt.x_true = x_true;
results_aug_adapt.w_true = w_true;
results_aug_adapt.x_hat = x_hat;
results_aug_adapt.w_hat = w_hat;
results_aug_adapt.eps_k_hist = eps_k_hist;

save('sim_results_aug_adapt.mat', '-struct', 'results_aug_adapt');
disp('Simulation complete. Adaptive results saved to sim_results_aug_adapt.mat');

% % 6. Quick Verification Plot
% figure('Name', 'Adaptive Filter Performance', 'Position', [100, 100, 800, 600]);

% subplot(2,1,1);
% plot(t, x_true(1,:), 'k', 'LineWidth', 1.5); hold on;
% plot(t, x_hat(1,:), 'g--', 'LineWidth', 1.5);
% ylabel('x_1'); title('State 1: True vs Adaptive Estimate');
% legend('True', 'Adaptive Augmented', 'Location', 'best');
% grid on;

% subplot(2,1,2);
% semilogy(t, eps_k_hist, 'b'); hold on;
% yline(tau_chi2, 'r--', 'Threshold (\tau)', 'LineWidth', 1.5);
% ylabel('NIS (\epsilon_k)'); xlabel('Time (s)');
% title('Normalized Innovation Squared (Log Scale)');
% grid on;