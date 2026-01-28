
x0=10;
b0=[0.3,-0.1,0.0]';
mu=0.5;

C = [1 1 0;
     1 0 1;
     1 0 0];

y = x0 + b0 + mu*randn(3,1);

theta_hat = (C.'*C)\(C.'*y);

x_hat  = theta_hat(1);
b_hat  = [theta_hat(2); theta_hat(3); 0];

