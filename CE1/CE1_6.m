clear all
close all
clc
%% Data preparation
N = 2000;
Ts = 0.25;
u = 0.7 * sign(rand(N, 1) - 0.5); % Binary signal between -0.7 and 0.7
T_sim = N * Ts;
t = (0:Ts:(N-1)*Ts)';

% Create struct simin for impulse response
simin.signals.values = u;
simin.time = t;

% Simulate impulse response
out_impulse = sim("CE1.slx");

% Extract impulse response as a vector y
y = out_impulse.simout.Data;

% Already prepare the segment length
segment_length = N / 8;
%% Whole data, no windowing
% 
R_yu_full = intcor(y, u);   
R_uu_full = intcor(u, u);

mid = floor(N/2) + 1;

% Only take the positive lag
R_yu = R_yu_full(mid:end-1);
R_uu = R_uu_full(mid:end-1);

% FFT without windowing
Y_whole = fft(R_yu);
U_whole = fft(R_uu);
G_whole = Y_whole ./ U_whole;

%% Windowing only

R_yu_full = intcor(y, u);   
R_uu_full = intcor(u, u);

mid = floor(N/2) + 1;

% Only take the positive lag
R_yu = R_yu_full(mid:mid+segment_length-1);
R_uu = R_uu_full(mid:mid+segment_length-1);

% Apply Hann window
h = hann(500);
h_shifted = h(251:end)';  % use second half 
R_yu_win = R_yu .* h_shifted;
R_uu_win = R_uu .* h_shifted;

% FFT
Y_win_whole = fft(R_yu_win);
U_win_whole = fft(R_uu_win);
G_win_whole = Y_win_whole ./ U_win_whole;
%% Averaging + Windowing

y_segments = cell(8,1);
u_segments = cell(8,1);  % Preallocate cell array to hold the segments

for i = 1:8
    idx_start = (i-1)*segment_length + 1;
    idx_end = i*segment_length;
    u_segments{i} = u(idx_start:idx_end);
    y_segments{i} = y(idx_start:idx_end);
end

u_fft = cell(8,1);
y_fft = cell(8,1);

num_lags = 500;
h = hann(num_lags);
mid = floor(num_lags / 2);
h_shifted = h(mid + 1:end)';

U_avg = zeros(1, segment_length);
Y_avg = zeros(1, segment_length);

for i = 1:8
    y_seg = y_segments{i};
    u_seg = u_segments{i};

    R_yu = xcorr(y_seg, u_seg, 'unbiased');
    R_uu = xcorr(u_seg, u_seg, 'unbiased');

    R_yu_tr = R_yu(mid : end)';
    R_uu_tr = R_uu(mid : end)';

    y_fft = fft(R_yu_tr .* h_shifted);
    u_fft = fft(R_uu_tr .* h_shifted);

    % Accumulate directly
    Y_avg = Y_avg + y_fft;
    U_avg = U_avg + u_fft;
end

Y_avg = Y_avg / 8;
U_avg = U_avg / 8;

% Final frequency response
G_avg = Y_avg ./ U_avg;

%% Plot

ws = 2*pi/(Ts);

segment_length = 1000;
freq_vect = (0:segment_length-1)' * (ws / segment_length);

% True system
G_true = tf([1.2], [1 2 1.35, 1.2]);

% Whole system
figure;
set(gcf, 'Position', [100, 100, 800, 400]);
sys_whole = frd(G_whole, freq_vect);
bode(sys_whole, 'b--');
hold on;
bode(G_true, 'r--');
xlim([0.1 30]);
ylim([-1440 0]);
title("Whole data");
legend('Whole data, no window', "True frequency response");

segment_length = 250;
freq_vect = (0:segment_length-1)' * (ws / segment_length);

% Windowing only
figure;
set(gcf, 'Position', [100, 100, 800, 400]);
sys_win_whole = frd(G_win_whole, freq_vect);
bode(sys_win_whole, 'b--');
hold on;
bode(G_true, 'r--');
xlim([0.1 30]);
ylim([-1440 0]);
title("Windowing only");
legend('Whole data, with window', "True frequency response");


% Windowing + averaging
figure;
set(gcf, 'Position', [100, 100, 800, 400]);
sys_win_avg = frd(G_avg, freq_vect);
bode(sys_win_avg, 'b--');
hold on;
bode(G_true, 'r--');
xlim([0.1 30]);
ylim([-1440 0]);
title("Windowing + averaging");
legend('Avg data, with window', "True frequency response");





