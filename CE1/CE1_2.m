clc
clear all
close all

% Define parameters
n = 6;  % Shift register length
p = 4;  % Number of periods

% Run the PRBS function
u = prbs(n, p);

% Compute autocorrelation using intcor function
[R, h] = intcor(u, u);

%%
% Plot the generated PRBS signal
figure;
stem(u, 'filled');
grid on;
title(['PRBS Signal with n=', num2str(n), ' and p=', num2str(p)]);
xlabel('Sample Index');
xlim([0 63]);
ylabel('Amplitude');

%%
% Plot the autocorrelation function
figure;
stem(h, R, 'filled');
grid on;
title('Autocorrelation Function of PRBS Signal');
xlabel('Lag (h)');
xlim([-150 150])
ylabel('R(h)');