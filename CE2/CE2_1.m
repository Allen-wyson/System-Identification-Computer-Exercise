clc; clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.1 Identification of a DC servomotor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.1.1 FIR model
clc; clear; close all;

path_ = "CE2_data.m";
[y, u, Ts] = GetExperimentData(path_);

N = size(u, 1); % input size
m = 200; % number of parameter, dimension of theta

tpc = u; % first column of toeplitz matrix
tpr = zeros(m, 1); % first row of toeplitz matrix
tpr(1,1) = u(1,1); % force 1st elt of input column to match first elt of 1st row

phi = toeplitz(tpc, tpr); % phi is a toeplitz matrix

theta_hat = phi \ y;    % estimated vector of parameters
y_hat = phi * theta_hat;    % predicted output
J = (norm(y - y_hat))^2 % loss function is the sum of squares of the prediction error

figure(1);
hold on;
plot(y);
plot(y_hat);
legend({'Measured output $y$','Predicted output $\hat{y}$'},'Interpreter','latex');
title('Output of the identified model');
grid;
xlabel('Time [ms]');

var_hat = 1 / (N - m) * J; % noise variance estimate
covar_hat = var_hat * inv(phi' * phi);
sigma = sqrt(diag(covar_hat));  % standard deviation of the parameter estimates

figure(2);
errorbar(theta_hat, 2 * sigma);
grid on;
title('Impulse response of the system with 95% confidence');
xlabel('Time [ms]')
ylabel('Amplitude')

%% 2.1.2 ARX to be continued
clc; clear; close all;

path_ = "CE2_data.m";
[y, u, Ts] = GetExperimentData(path_);

N = size(u,1); %number of input

function [Big_Phi] = create_BigPhi_matrix(u, y, N)
        Big_Phi = zeros(N,4);
        Big_Phi(2,:) = [-y(1), 0, u(1), 0];
        for k=3:N
            Big_Phi(k,:)= [-y(k-1),-y(k-2),u(k-1),u(k-2)];
        end
    end

BigPhi = create_BigPhi_matrix(u, y, N);

%Question 1
theta_ARX = inv(BigPhi' * BigPhi) * BigPhi' * y % estimated vector of parameters

%Question 2
y_ARX = BigPhi * theta_ARX;
J_ARX = norm(y - y_ARX) ^ 2 % loss function is the sum of squares of the prediction error

figure;
hold on;
plot(y, 'b');
plot(y_ARX, 'r--');
legend("True", "ARX Predicition")
title("Comparison");


%Question 3
a1 = theta_ARX(1);
a2 = theta_ARX(2);
b1 = theta_ARX(3);
b2 = theta_ARX(4);

f_sampling = 1e3;

q = tf('q', 1/f_sampling);
sys = (b1*(1/q) + b2*(1/q^2)) / ...
        (1 + a1*(1/q) + a2*(1/q^2));

% Simulate transfer function model
t = 0:(1/f_sampling):(N-1)*(1/f_sampling);
ym = lsim(sys, u, t);

J_tf = norm(ym - y)^2

figure;
hold on;
plot(y, 'b');
plot(ym, 'r--');
legend("True", "ym")
title("Comparison");

% Question 4 - Instrumental Variable (IV) Method

% Step 1: Create the IV matrix (using y_ARX as instrumental variables)
Big_Phi_IV = zeros(N,4);
Big_Phi_IV(2,:) = [-y_ARX(1), 0, u(1), 0];
for k=3:N
    Big_Phi_IV(k,:) = [-y_ARX(k-1), -y_ARX(k-2), u(k-1), u(k-2)];
end

% Step 2: Estimate parameters with IV
theta_IV = inv(Big_Phi_IV' * BigPhi) * Big_Phi_IV' * y  % IV estimator formula

% Step 3: Compute predicted output with IV model
y_IV = BigPhi * theta_IV;
J_IV = norm(y - y_IV)^2  % Loss function for IV

% Step 4: Compare with ARX results
fprintf('ARX Loss: %.4f | IV Loss: %.4f\n', J_ARX, J_IV);

% Step 4: Plot ARX vs IV in subplots
figure;
subplot(1,2,1); % Left subplot: ARX
hold on;
plot(y, 'b', 'LineWidth', 1.5);
plot(y_ARX, 'r--', 'LineWidth', 1.5);
legend("True", "ARX Prediction");
title("ARX Model");
xlabel("Time (samples)"); ylabel("Output");
grid on;

subplot(1,2,2); % Right subplot: IV
hold on;
plot(y, 'b', 'LineWidth', 1.5);
plot(y_IV, 'g:', 'LineWidth', 1.5);
legend("True", "IV Prediction");
title("IV Model");
xlabel("Time (samples)"); 
grid on;

% Step 5: Simulate IV model (like Q3 for ARX)
sys_IV = (theta_IV(3)*(1/q) + theta_IV(4)*(1/q^2)) / ...
          (1 + theta_IV(1)*(1/q) + theta_IV(2)*(1/q^2));
ym_IV = lsim(sys_IV, u, t);

% Compare errors
J_tf_IV = norm(ym_IV - y)^2;
fprintf('ARX Sim Error: %.4f | IV Sim Error: %.4f\n', J_tf, J_tf_IV);

%% 2.1.3 state-space
clc; clear; close all;

path_ = "CE2_data.m";
[y, u, Ts] = GetExperimentData(path_);

r = 10; % maximal order of the system
N = size(u, 1); % input size

%Question 1

% Construct U and Y
Y = zeros(r, N);
U = zeros(r, N);
for i=1:N
	for j=1:r
		k = i+j-1;
		if k > N
			Y(j,i) = 0;
			U(j,i) = 0;
		else
			Y(j,i) = y(k);
			U(j,i) = u(k);
		end
	end
end

U_orthogonal = eye(N) - U' * inv(U*U') * U;
Q = Y * U_orthogonal;

%Question 2
singular_value = svd(Q);
plot(singular_value);
[UU, S, V] = svd(Q);

% from the plot we guess n = 3
n = 3;
Or = UU(:, 1:n);

%Question 3
C = Or(1, :);
A_hat = pinv(Or(1:(r-1),:)) * Or(2:r,:);

% compute uf and then B estimate
q = tf('z');
F = C * inv(q * eye(n) - A_hat);
uf = zeros(N, n);
for i = 1 : n
    uf(:, i) = lsim(F(i), u);
end

% assume D = 0
D = 0;
B_hat = pinv(uf) * y

%Question 4
sys = ss(A_hat, B_hat, C, D, Ts); %add 1e-3 as sapling time to have a discrete system
[y_sys, t_sys, x_sys] = lsim(sys, u);

fprintf('Loss function: %d\n', norm(y_sys - y)^2);

figure
plot(t_sys, y_sys); % simulated output ym
hold on
t = linspace(0, (size(y,1)-1)*Ts, size(y,1));
plot(t, y, 'r') %real output y
legend("simulated output ym", "real output y")
grid

figure
plot(t_sys, x_sys);
grid;
legend("x_1", "x_2", "x_3") % assumed order n of the system = 3
