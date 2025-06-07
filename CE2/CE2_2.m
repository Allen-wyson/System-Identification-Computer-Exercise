%% 2.2.1
% Question 1

clc; clear; close all;
load data_position1.mat

% Compute the derivative of the output signal
y_derivative = lsim(1 - tf('z', Ts)^-1, y);

% Create detrended IDDATA object
Z = detrend(iddata(y_derivative, u, Ts, "Period", 8191));
nk = 1; % Initial delay
loss_values = zeros(1, 15);
models = cell(1, 15);

% Identify 15 ARX models with orders from 1 to 15
for delta = 1:15

disp(['Order = ', num2str(delta)]);
sys = arx(Z, [delta delta nk]);
loss_values(delta) = sys.EstimationInfo.LossFcn;
models{delta} = sys; % Store each model for later use

end

%% Plot the loss function evolution

x_axis = 1:15;
figure;
plot(x_axis, loss_values, '-o');
xlabel('Model order \delta');
ylabel('Loss function value');
title('Loss Function vs. Model Order');
grid on;

%% Plot the Bode diagram of the best model

% Select the model with the minimum loss value
best_order = 8;
best_model = models{best_order};

% Plot the Bode diagram
figure;
bode(best_model);
title(['Bode Diagram of Best ARX Model (Order = ' num2str(best_order) ')']);
grid on;

%% Question 2

close all;
delta_selected = [5, 6, 7, 8, 9, 10];
figure();

for i = 1:length(delta_selected)

sys = armax(Z, [delta_selected(i) delta_selected(i) delta_selected(i) 1]);
subplot(2, 3, i);
h = iopzplot(sys);
hold on;
showConfidence(h, 2);
title(['Order = ' num2str(delta_selected(i))]);

end

%% Question 3
% Fix model order (based on previous steps)
na = 8; % You can use the estimated order from the loss function curve
nb = 8;
nc = 8;
fprintf('\nEstimating nk by checking B coefficients:\n');

% Fit ARMAX model
model = armax(Z, [na nb nc 0]);
% Get numerator coefficients and their std deviations
b = model.b;
db = model.db;
lower = b - 2*db;
upper = b + 2*db;
idx = (lower <= 0) & (upper >= 0) % 1×8 logical vector
%nk = find(idx, 1); % the indices k where 0∈[b_k−2db_k, b_k+2db_k]
nk = 1; %TBD by looking at idx first values which are 1

% Display and analyze
fprintf('nk = %d | B = %s | std = %s\n', nk, mat2str(b, 3), mat2str(db, 3));

% Optional: plot coefficients vs. std deviation
figure;
errorbar(1:length(b), b, db, 'o');
title(['B Coefficients with Std Dev for nk = ' num2str(nk)]);
xlabel('Coefficient Index');
ylabel('Value');
grid on;

%% Question 4 – Comparison with struc, arxstruc, selstruc

% Split the data into estimation and validation sets
N = length(Z.y);
Z_est = Z(1:round(2*N/3)); % 2/3 for estimation
Z_val = Z(round(2*N/3)+1:end); % 1/3 for validation

% Define model structures to test: na = 1:15, nb = 1:15, nk = 1 (fixed)
orders = struc(1:15, 1:15, 1);
% Evaluate ARX models on the estimation and validation sets
V = arxstruc(Z_est, Z_val, orders);
% Automatically select the best model using FPE (default)
[best_order, best_index] = selstruc(V, 0);
% Display the selected model structure
fprintf('\nBest ARX structure according to selstruc:\n');
fprintf('na = %d, nb = %d, nk = %d\n', best_order);
% Estimate ARX model with the best selected structure
best_model = arx(Z, best_order);
% Print the loss function of this model
fprintf('LossFcn for selected model = %.4f\n', best_model.EstimationInfo.LossFcn);


%% 2.2.2 – Parametric Identification
% Reload data if needed
load data_position1.mat

% Compute the derivative of the output signal
y_derivative = lsim(1 - tf('z', Ts)^-1, y);
y_derivative = detrend(y_derivative, 'constant');
u = detrend(u, 'constant');

% Create detrended IDDATA object
Z = detrend(iddata(y_derivative, u, Ts, "Period", 8191));

% Split data into estimation and validation sets
N = length(Z.y);
Z_est = Z(1:round(2*N/3)); % Use 2/3 of data for estimation
Z_val = Z(round(2*N/3)+1:end); % Use 1/3 of data for validation

% Define model orders based on previous identification
na = 8; % Number of past outputs in the model A(q^-1)
nb = 8; % Number of past inputs in the model B(q^-1)
nk = 1; % Input delay (number of samples before input affects output)
nc = na; % Order of C(q^-1) in ARMAX/BJ (noise model)
nd = na; % Order of D(q^-1) in BJ (noise model)
nf = nb; % Order of F(q^-1) in OE/BJ (denominator of dynamic model)

% Estimate models using different structures
models = struct();
models.arx = arx(Z_est, [na nb nk]); % ARX: A(q^-1)y(t) = B(q^-1)u(t-nk)
models.iv4 = iv4(Z_est, [na nb nk]); % IV4: Instrumental variable version of ARX
models.armax = armax(Z_est, [na nb nc nk]); % ARMAX: Adds C(q^-1) noise model
models.oe = oe(Z_est, [nb na nk]); % OE: Output-Error model with B(q^-1)/F(q^-1)
models.bj = bj(Z_est, [nb na na na nk]); % BJ: Box-Jenkins model with full noise model
models.n4sid = n4sid(Z_est, na); % N4SID: State-space model with total order

%% 2.2.3.1 Model Validation
modelNames = fieldnames(models);    % {'arx','iv4','armax','oe','bj','n4sid'}
nModels    = numel(modelNames);

figure('Name','Model Validation','Units','normalized','Position',[.1 .1 .8 .7]);

for k = 1:nModels
    sys = models.(modelNames{k});   % extract the k-th idmodel
    subplot(3,2,k);
    
    % Plot measured vs. simulated output & show fit (%)
    compare(Z_val, sys);
    % k
    % modelNames
    % Title with model name and fit %
    % (you can get the fit value if you need it numerically:
    % [~, fit, ~] = compare(Z_val,sys); )
    title(sprintf('%s Model', modelNames{k}));
end

% improve spacing
% tightfig;    % if you have tightfig, or use:
set(gcf,'Renderer','painters');

%% 2.2.3.2 freq response
G_spa = spa(Z_est);
figure;
hold on;
bode(models.arx)
bode(models.iv4)
bode(models.armax)
bode(models.oe)
bode(models.bj)
bode(models.n4sid)
bode(G_spa)
hold off;
legend('arx', 'iv4', 'armax', 'oe', 'bj', 'n4sid', 'True model')

%% 2.2.3.3 whiteness test
figure;
resid(Z_val, models.arx); title('ARX');
figure;
resid(Z_val, models.iv4); title('Intrumental Variables');
figure;
resid(Z_val, models.armax); title('ARMAX');
figure;
resid(Z_val, models.oe); title('Output Error');
figure;
resid(Z_val, models.bj); title('Box Jenkins');
figure;
resid(Z_val, models.n4sid); title('State-Space');

%% 2.2.3.4 — Add integrator and re‐validate on second half of position data

% 1) Reload raw position data (y) and input (u)
load data_position1.mat   % brings in u, y, Ts

% 2) Build iddata for the SECOND HALF of the raw position signal
Npos = length(y);
idx2 = round(Npos/2)+1 : Npos;
y_second = y(idx2);
u_second = u(idx2);
posValData = iddata(y_second, u_second, Ts);

% 3) Build the discrete‐time integrator 1/(1 - z^{-1}) as an idtf
integrator = idtf(1, [1 -1], Ts);

% 4) Loop over each derivative‐domain model, append integrator, and compare
figure('Name','Position‐Domain Validation','Units','normalized','Position',[.1 .1 .8 .7]);
for k = 1:nModels
    sysDer = models.(modelNames{k});       % derivative‐domain model
    sysPos = sysDer * integrator;          % now maps u → position

    subplot(3,2,k);
    compare(posValData, sysPos);
    title(sprintf('%s Model', modelNames{k}));
end
set(gcf,'Renderer','painters');

% 5) Nonparametric FRF of u→position using SPA on posValData
G_pos_spa = spa(posValData);

% 6) Overlay Bode plots: six integrated models vs. nonparametric FRF
figure('Name','Bode: Parametric vs. Nonparametric (Position)','Units','normalized','Position',[.1 .1 .8 .7]);
hold on;
for k = 1:nModels
    sysDer = models.(modelNames{k});
    sysPos = sysDer * integrator;
    bode(sysPos);
end
bode(G_pos_spa);
hold off;
legend('arx', 'iv4', 'armax', 'oe', 'bj', 'n4sid', 'True model')
