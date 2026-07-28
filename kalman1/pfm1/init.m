%% 1. Define Continuous-Time System Matrices
Ac = [1, 2; -3, -4];
Buc = [1; 0];
Bwc = [0.4; 0.1];
C = [1, 2; 1, 0; 1, 0];
Du = [0; 0; 0];
Dw = [1; 2; 3]; 

%% 2. Discretization using zero-order hold
tau = 1e-4;
Ad=eye(size(Ac,1))+tau*Ac;
Bud=tau*Buc;
Bwd=tau*Bwc;

n_x = size(Ad, 1);       % Number of physical states
n_w = size(Bwd, 2);      % Number of disturbances
n_y = size(C, 1);        % Number of measurements

N = 1000;              % Number of simulation steps
t = (0:N-1) * tau;     % Time vector
u = sin(2*pi*5*t);     % Example control input

% 3. Noise Covariances
Q_x = diag([1e-6, 1e-6]); 
Q_w = 1e-4;               
R   = diag([1e-4, 1e-4, 1e-4]); 

% augmented kalman favoring disturbance
w_true = zeros(n_w, N);
w_true(1, 201:400) = 2.0;                                     
w_true(1, 401:700) = 2.0 + 1.5 * sin(2*pi*2*t(401:700));      
w_true(1, 701:N)   = linspace(2.0, -1.0, N-700);              

% % dual psm favoring disturbance
% w_true = zeros(n_w, N);
% w_true(1, 201:500) = 5 * sign(sin(2*pi*50*t(201:500))); 
% w_true(1, 600) = 15; 
% w_true(1, 700) = -20;
% w_true(1, 800) = 25;
% w_true(1, 900) = -15;