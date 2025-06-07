clear all;
Ts = 0.25;
T_sim = 100;
N = T_sim / Ts;

%% **Exercise 1: Step Response**
% Define simulation time
T_end = 100; % Total simulation time in seconds
t = (0:Ts:(N-1)*Ts)'; % Generates exactly 400 elements

% Define step input
u_step = zeros(size(t)); % Initialize with zeros
u_step(t >= 1) = 1; % Step occurs at t = 1s

% Create struct simin for step response
simin.signals.values = u_step;
simin.time = t;

% Simulate step response
out_step = sim("CE1.slx");

% Plot step response
figure;
plot(out_step.simout.Time, out_step.simout.Data, 'b', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Response');
title('Step Response');
legend('System Output');

%% **Exercise 1: Impulse Response**
% Define impulse input (delta function)
u_impulse = zeros(size(t)); % Initialize with zeros
u_impulse(t == 1) = 1; % Impulse at t = 1

% Create struct simin for impulse response
simin.signals.values = u_impulse;
simin.time = t;

% Simulate impulse response
out_impulse = sim("CE1.slx");

% Extract impulse response as a vector y
y = out_impulse.simout.Data;
size(y)

% Plot impulse response
figure;
plot(out_impulse.simout.Time, y, 'r', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Response');
title('Impulse Response');
legend('System Output');