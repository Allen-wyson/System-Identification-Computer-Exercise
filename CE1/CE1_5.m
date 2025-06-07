clear all;
close all;
clc;
%% Part 1
Ts = 0.25;
u = prbs(8, 8) *0.7;
N = size(u,1);
T_sim = N * Ts;
t = (0:Ts:(N-1)*Ts)';

% Create struct simin for impulse response
simin.signals.values = u;
simin.time = t;
%% Part 2

% Simulate impulse response
out_impulse = sim("CE1.slx");

% Extract impulse response as a vector y
y = out_impulse.simout.Data;

%% Part 3

y_segments = cell(8,1);
u_segments = cell(8,1);  % Preallocate cell array to hold the segments
segment_length = N / 8;

for i = 1:8
    idx_start = (i-1)*segment_length + 1;
    idx_end = i*segment_length;
    u_segments{i} = u(idx_start:idx_end);
    y_segments{i} = y(idx_start:idx_end);
end

u_fft = cell(8,1);
y_fft = cell(8,1);

for i = 1:8
    u_fft{i} = fft(u_segments{i});
    y_fft{i} = fft(y_segments{i});
end




% Initialize average vector with zeros of the same size as one segment
U_avg = zeros(segment_length, 1);
Y_avg = zeros(segment_length, 1);

% Sum over all segments
for i = 2:8
    U_avg = U_avg + u_fft{i};
    Y_avg = Y_avg + y_fft{i};
end

% Divide by 8 to get the average
U_avg = U_avg / 7;
Y_avg = Y_avg / 7;

G_avg = Y_avg ./ U_avg;

%%
ws = 2*pi/(Ts);
freq_vect = zeros(segment_length, 1);

for i=1:(segment_length-1)
    freq_vect(i+1) = ws * i / segment_length;
end

sys = frd(G_avg, freq_vect);

bode(sys);

hold on

G_true = tf([1.2], [1 2 1.35, 1.2])
bode(G_true, 'r--');

hold off
legend('Identified', 'True');

