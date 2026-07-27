% 1. Define Continuous-Time System Matrices
Ac = [1, 2; -3, -4];
Buc = [1; 0];
Bwc = [0.1; 0.1];
C = [1, 2; 1, 0; 1, 0];
Du = [0; 0; 0];
Dw = [1; 2; 3]; 

% 2. Discretization using zero-order hold
tau = 1e-3;
Ad=eye(size(Ac,1))+tau*Ac;
Bud=tau*Buc;
Bwd=tau*Bwc;

n_x = size(Ad, 1);       % Number of physical states
n_w = size(Bwd, 2);      % Number of disturbances
n_y = size(C, 1);        % Number of measurements

