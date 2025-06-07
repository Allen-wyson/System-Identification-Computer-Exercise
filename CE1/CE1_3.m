clc; clear all; close all;

Ts = 0.25;
T_sim = 100;
N = T_sim / Ts + 1;

u = (rand(N, 1) - 0.5) * 1.4;    % random signal between -0.7 and 0.7
%u = 0.7 * rand(N, 1);


t = 0 : Ts: (N - 1) * Ts; % time vector

simin.signals.values = u;
simin.time = t';

out = sim("CE1.slx");
y = out.simout.Data;

L = 100; % assumed length of the impulse response
U = toeplitz(u, [u(1); zeros(L - 1, 1)]);
Theta_K = inv(U'*U) * U' * y(1:size(U, 1));

lambda = 0.1; % Regularization parameter (adjust as needed)
I = eye(N); % Identity matrix of size K

% Compute regularized impulse response
U = toeplitz(u, [u(1); zeros(N - 1, 1)]);
Theta_K_reg = inv(U' * U + lambda * I) * U' * y(1:size(U, 1));

% Define the continuous-time transfer function
Gs = tf(1.2, [1 2 1.35 1.2]);

% Discretize with zero-order hold (ZOH)
Gz = c2d(Gs, Ts, 'zoh');
true_impulse_response = impulse(Gz, t)*Ts;

%% ERROR NORM 2

% Calculation of the error for the finite impulse response
error_finite = Theta_K - true_impulse_response(1:L);
norm_error_finite = norm(error_finite, 2);  % 2-norm of the error for h_finite

% Calculation of the error for the regularized impulse response
error_regularized = Theta_K_reg(1:L) - true_impulse_response(1:L);
norm_error_regularized = norm(error_regularized, 2);  % 2-norm of the error for h_regularized

% Displaying the results
disp(['2-norm of the error for the finite impulse response: ', num2str(norm_error_finite)]);
disp(['2-norm of the error for the regularized impulse response: ', num2str(norm_error_regularized)]);



%% Plot
figure(1)
hold on;
plot(Theta_K,'b');
plot(Theta_K_reg(1: L),'r');
plot(true_impulse_response(1:L), 'k');
xlabel('Sample Index k');
ylabel('Impulse Response g(k)');
legend('Finite Impulse Response', 'Regularized Impulse Response', 'True Impulse Response')
title('Impulse Response');
grid on;
