% plot_comparisons.m
clear; clc; close all;

% 1. Load the saved datasets
% These will load as structures containing t, x_true, w_true, x_hat, w_hat
data_uio = load('sim_results_uio.mat');
data_aug = load('sim_results_aug.mat');

% Verify that both time vectors are identical (optional safety check)
if ~isequal(data_uio.t, data_aug.t)
    warning('Time vectors do not match. Ensure both scripts ran with the same parameters.');
end

t = data_uio.t;

% 2. Plotting
figure('Name', 'Filter Comparison: Dual-PFM vs Augmented', 'Position', [100, 100, 800, 800]);

% --- State 1 ---
subplot(3,1,1);
plot(t, data_uio.x_true(1,:), 'k', 'LineWidth', 1.5); hold on;
plot(t, data_uio.x_hat(1,:), 'r--', 'LineWidth', 1.5);
plot(t, data_aug.x_hat(1,:), 'b--', 'LineWidth', 1.5);
ylabel('x_1'); 
title('State 1 Estimation');
legend('True State', 'Dual-PFM Estimate', 'Augmented Estimate', 'Location', 'best');
grid on;

% --- State 2 ---
subplot(3,1,2);
plot(t, data_uio.x_true(2,:), 'k', 'LineWidth', 1.5); hold on;
plot(t, data_uio.x_hat(2,:), 'r--', 'LineWidth', 1.5);
plot(t, data_aug.x_hat(2,:), 'b--', 'LineWidth', 1.5);
ylabel('x_2'); 
title('State 2 Estimation');
legend('True State', 'Dual-PFM Estimate', 'Augmented Estimate', 'Location', 'best');
grid on;

% --- Disturbance ---
subplot(3,1,3);
plot(t, data_uio.w_true(1,:), 'k', 'LineWidth', 1.5); hold on;
plot(t, data_uio.w_hat(1,:), 'r--', 'LineWidth', 1.5);
plot(t, data_aug.w_hat(1,:), 'b--', 'LineWidth', 1.5);
ylabel('w'); xlabel('Time (s)'); 
title('Disturbance Estimation');
legend('True Disturbance', 'Dual-PFM Estimate', 'Augmented Estimate', 'Location', 'best');
grid on;

% 3. Calculate and Display Mean Squared Error (MSE)
fprintf('--- Mean Squared Error (MSE) Comparison ---\n');

% Dual-PFM Errors
mse_uio_x1 = mean((data_uio.x_true(1,:) - data_uio.x_hat(1,:)).^2);
mse_uio_x2 = mean((data_uio.x_true(2,:) - data_uio.x_hat(2,:)).^2);
mse_uio_w  = mean((data_uio.w_true(1,:) - data_uio.w_hat(1,:)).^2);

% Augmented Errors
mse_aug_x1 = mean((data_aug.x_true(1,:) - data_aug.x_hat(1,:)).^2);
mse_aug_x2 = mean((data_aug.x_true(2,:) - data_aug.x_hat(2,:)).^2);
mse_aug_w  = mean((data_aug.w_true(1,:) - data_aug.w_hat(1,:)).^2);

fprintf('State 1 (x1):\n');
fprintf('  Dual-PFM MSE: %e\n', mse_uio_x1);
fprintf('  Augmented MSE: %e\n', mse_aug_x1);

fprintf('State 2 (x2):\n');
fprintf('  Dual-PFM MSE: %e\n', mse_uio_x2);
fprintf('  Augmented MSE: %e\n', mse_aug_x2);

fprintf('Disturbance (w):\n');
fprintf('  Dual-PFM MSE: %e\n', mse_uio_w);
fprintf('  Augmented MSE: %e\n', mse_aug_w);