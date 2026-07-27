decay_factor = 0.9999; % Keeps the disturbance mode strictly detectable

A_aug = [Ad, Bwd; 
         zeros(n_w, n_x), decay_factor * eye(n_w)]; 
C_aug = [C, Dw];


n_z = size(A_aug, 1);
n_y = size(C_aug, 1);

%% 2. Define Noise Weights
W_w = 1.0;   % Process noise variance (disturbance random walk)
W_v = 0.05;  % Sensor noise variance

% Map the noise inputs to the augmented states and outputs
B_noise = [zeros(n_x, n_w); 
           sqrt(W_w) * eye(n_w)];
D_noise = sqrt(W_v) * eye(n_y);

%% 3. Formulate the H2 Optimization in CVX
% We seek to minimize Trace(Z), where Z bounds the error covariance.
% The error dynamics are: e(k+1) = (A_aug - L*C_aug)e(k) + B_noise*w_v - L*D_noise*v
cvx_begin sdp quiet
    % Variables
    variable P(n_z, n_z) symmetric
    variable Y(n_z, n_y)
    variable Z(n_z, n_z) symmetric
    
    % Objective: Minimize the H2 performance bound
    minimize( trace(Z) )
    
    subject to
        % Constraint 1: Covariance bound Z > P^-1 
        % (Schur complement of Z - P^-1 > 0)
        [ Z,        eye(n_z);
          eye(n_z), P ] >= 1e-6 * eye(2*n_z);
        
        % Constraint 2: Discrete-time H2 Stability and Performance LMI
        % This is the Schur complement of the discrete Lyapunov inequality 
        % incorporating the disturbance matrices.
        [ P,             P*A_aug - Y*C_aug,  P*B_noise,       -Y*D_noise;
          A_aug'*P - C_aug'*Y', P,           zeros(n_z, n_w), zeros(n_z, n_y);
          B_noise'*P,    zeros(n_w, n_z),    eye(n_w),        zeros(n_w, n_y);
          -D_noise'*Y',  zeros(n_y, n_z),    zeros(n_y, n_w), eye(n_y) ] >= 1e-6 * eye(2*n_z + n_w + n_y);
          
        % Optional: Regional Pole Placement 
        % Uncomment to enforce a maximum decay rate (e.g., r = 0.98) 
        % to prevent the H2 solver from picking excessively slow poles.
        % r = 0.98;
        % [ r * P,               P * A_aug - Y * C_aug;
        %   A_aug'* P - C_aug'* Y', r * P ] >= 1e-6 * eye(2*n_z);
          
cvx_end

%% 4. Recover Gain and Verify
if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
    % Recover the true observer gain
    L_luen = P \ Y; 
    
    disp('Success: H2-Optimal observer gain L_luen designed via CVX.');
    fprintf('Optimal H2 Upper Bound (Trace(Z)): %.4f\n', trace(Z));
    
    % Verify error dynamics
    A_err = A_aug - L_luen * C_aug;
    poles_err = eig(A_err);
    fprintf('Max closed-loop error pole magnitude: %.4f\n', max(abs(poles_err)));
else
    disp('CVX Status:');
    disp(cvx_status);
    error('CVX could not find a feasible solution.');
end