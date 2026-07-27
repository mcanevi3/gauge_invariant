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

% 3. H2 Filter Tuning Weights
% These are no longer "variances", just tuning knobs.
% Increase W to track the bump better; increase V_noise to trust sensors less.
W = 100.0;       
V_noise = 0.05; 

Q = Bwd * W * Bwd'; 
R = (Dw * W * Dw') + (V_noise * eye(3)); 
N = Bwd * W * Dw';  

% 4. Initialization
steps = 1500; % Extended to see the sine wave clearly
x_true = [0; 0];
x_est  = [0; 0];
P      = eye(2);
u      = zeros(1, steps); % Set control input to zero to just see disturbance response

x_true_hist = zeros(2, steps);
x_est_hist  = zeros(2, steps);

% Define the deterministic disturbance (e.g., a 1 Hz sine wave road profile)
t = (0:steps-1) * Ts;
w_deterministic = 2 * sin(2 * pi * 5.0 * t); 

% 5. Optimal Observer Loop (H2 Filter)
for k = 1:steps
    % --- A. SYSTEM SIMULATION ---
    w_k = w_deterministic(k);          % Deterministic physical disturbance
    v_k = sqrt(V_noise) * randn(3, 1); % Sensors still have random static
    
    y_k = C * x_true + Dw * w_k + v_k;
    
    x_true_hist(:, k) = x_true;
    x_est_hist(:, k)  = x_est; 
    
    % --- B. H2 GAIN CALCULATION ---
    S = C * P * C' + R;
    L = (Ad * P * C' + N) / S; 
    
    % --- C. PREDICT NEXT STATE ---
    x_est_next = Ad * x_est + Bud * u(k) + L * (y_k - C * x_est);
    P_next = Ad * P * Ad' + Q - L * S * L';
    
    % --- D. ADVANCE TRUE SYSTEM ---
    x_true = Ad * x_true + Bud * u(k) + Bwd * w_k;
    
    x_est = x_est_next;
    P = P_next;
end

% 6. Plot the Results
figure(1);
subplot(2,1,1);cla;hold on;grid on;
xlabel('t');ylabel('x_1');title('Augmented Filter: State and Disturbance Estimation'); 
plot(t, x_true_hist(1,:), 'k', 'LineWidth', 3);
plot(t, x_est_hist(1,:), 'r', 'LineWidth', 2);
legend('True x_1', 'Estimated x_1'); 

subplot(2,1,2);cla;hold on;
xlabel('t');ylabel('x_2'); grid on;
plot(t, x_true_hist(2,:), 'k', 'LineWidth', 3);
plot(t, x_est_hist(2,:), 'r', 'LineWidth', 2);
legend('True x_2', 'Estimated x_2'); 
