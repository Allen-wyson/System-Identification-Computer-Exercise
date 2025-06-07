clear all;
close all;
clc;
%% 1.4.1 Generate input sequence, apply to the system

Ts = 0.25;
u = prbs(7, 2)*0.7;
N = size(u,1);
T_sim = N * Ts;
t = (0:Ts:(N-1)*Ts)';

% Create struct simin for impulse response
simin.signals.values = u;
simin.time = t;

% Simulate impulse response
out_impulse = sim("CE1.slx");

% Extract impulse response as a vector y
y = out_impulse.simout.Data;

% Plot impulse response
% figure;
% plot(out_impulse.simout.Time, y, 'r', 'LineWidth', 2);
% grid on;
% xlabel('Time (s)');
% ylabel('Response');
% title('Impulse Response');
% legend('System Output');


%% Correlation approach

Ruu = intcor(u, u);
Ryu = intcor(y, u);

% Resize the vectors
Ruu = Ruu(1:N);
Ryu = Ryu(1:N);

% Construct Toeplitz input matrix
row = [Ruu(1), zeros(1, N-1)];
col = Ruu(1:N);
RUU = toeplitz(col, row);

K = 80; % 20*4, estimation of the response time

% Troncation
RUUk = RUU(1:N, 1:K);

% Impulse response computation
imp_res = inv(RUUk'*RUUk) * RUUk' *Ryu';

% figure
% % Subplot 1: Non-regularized impulse response
% stem(imp_res, 'filled', 'b');
% xlabel('Sample Index k');
% ylabel('Impulse Response g(k)');
% title('Non-Regularized Impulse Response');
% grid on;



%% Correlation with matlab function and comparison

Ruu_matlab = xcorr(u,u, 'biased');
Ryu_matlab = xcorr(y,u, 'biased');

Ruu_matlab = Ruu_matlab(1:N);
Ryu_matlab = Ryu_matlab(1:N);

% Construct Toeplitz input matrix
row = [Ruu_matlab(1), zeros(1, N-1)];
col = Ruu_matlab(1:N);
RUU_matlab = toeplitz(col, row);

% Troncation
RUUk_matlab = RUU_matlab(1:N, 1:K);

% Impulse response computation
imp_res_matlab = inv(RUUk_matlab'*RUUk_matlab) * RUUk_matlab' *Ryu_matlab;

%% Plots

% Impulse
u = zeros(N,1);
u(1) = 1; % Impulse at first sample

% Create struct simin for impulse response
simin.signals.values = u;
simin.time = t;

% Simulate impulse response
out_impulse = sim("CE1_without_noise.slx");

% Extract impulse response as a vector y
y_imp = out_impulse.simout.Data;

% Compute 2-norm
g_true = y_imp(1:K);

% Compute 2-norm errors
error_intcor = norm(imp_res - g_true);
error_xcorr  = norm(imp_res_matlab - g_true);





figure
hold on
h1 = stem(imp_res, 'filled', 'b', 'DisplayName', 'From intcor()');
h2 = stem(g_true, 'filled', 'g', 'DisplayName', 'True Impulse Response');
h3 = stem(imp_res_matlab, 'filled', 'r', 'DisplayName', 'From xcorr()');
hold off

xlabel('k');
ylabel('Impulse Response g(k)');
title('Impulse Response Comparison');
legend('Location', 'best');
grid on;

% Print errors
fprintf('2-norm error (intcor): %.4f\n', error_intcor);
fprintf('2-norm error (xcorr):  %.4f\n', error_xcorr);







