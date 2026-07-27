
% 7. Plot the Results
figure(1);
subplot(3,1,1);cla;hold on;grid on;
xlabel('t');ylabel('x_1');title('State and Disturbance Estimation'); 
plot(t, x_true_hist(1,:), 'k', 'LineWidth', 3);
plot(t, x_kf_hist(1,:), 'r', 'LineWidth', 2);
plot(t, x_luen_hist(1,:), 'b', 'LineWidth', 2);
legend('x_1', 'Kalman','Luenberger'); 

subplot(3,1,2);cla;hold on;
xlabel('t');ylabel('x_2'); grid on;
plot(t, x_true_hist(2,:), 'k', 'LineWidth', 3);
plot(t, x_kf_hist(2,:), 'r', 'LineWidth', 2);
plot(t, x_luen_hist(2,:), 'b', 'LineWidth', 2);
legend('x_2', 'Kalman','Luenberger'); 

subplot(3,1,3);cla;hold on;grid on;
xlabel('t'); ylabel('w');
plot(t, x_true_hist(3,:), 'k', 'LineWidth', 3);
plot(t, x_kf_hist(3,:), 'r', 'LineWidth', 2);
plot(t, x_luen_hist(3,:), 'b', 'LineWidth', 2);
legend('w', 'Kalman','Luenberger'); 
