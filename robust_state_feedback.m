% Robust State-Feedback Controller Design for Two-Cart and Spring System
% System matrices
close all
clear all

A = [0.5780 0.8492 0.4220 0.1508;
     -0.6985 0.5780 0.6985 0.4220;
     0.4220 0.1508 0.5780 0.8492;
     0.6985 0.4220 -0.6985 0.5780];

B = [0.4610; 0.8492; 0.0390; 0.1508];

% Simulation parameters
t = 1; % Lag
n = 4; % Number of states
m = 1; % Number of inputs
L = n + t;
T = 100; % Number of time steps
d_max = 0.1; % Maximum disturbance magnitude
mag = 1;




% Generate persistently exciting input
u = zeros(m, T);
u = (2*mag).*rand(m,T)-mag;

% Simulate system with noise
x = zeros(n, T+1);
x(:, 1) = (2*mag).*rand(n,1)-mag;
D0 = zeros(n, T);

for k = 1:T
    d = d_max * (2*rand(n,1) - 1); % Random disturbance
    D0(:,k) = d;
    x(:,k+1) = A*x(:,k) + B*u(:,k) + d;
end

% Construct data matrices
X0 = x(:,1:T);
X1 = x(:,2:T+1);
U0 = u;

if rank([U0 ; X0]) == m+n
disp('data are sufficiently rich');
end

% Check if DD' < I*delta
if D0*D0' <= eye(n)*d_max^2
disp('delta bounded');
end
Theta = diag([0 0 1 0]);

% Solving the non-robust state feedback controller
cvx_begin sdp

    variable Y(T,n)
    variable S(n,n) symmetric

    [S-eye(n) ((X1-D0)*Y); ((X1-D0)*Y)' S] >= 0;
    S == X0*Y;
    
cvx_end

K = U0*Y/S;


if strcmpi(cvx_status, 'Solved')
    disp('Robust controller design successful');
else
    error('Robust controller design failed');
end 

% Simulate closed-loop system
x_cl = zeros(n, T+1);
U1 = zeros(m, T);
x_cl(:,1) = x(:,1); % Use the same initial condition
u_ctl = zeros(n, T+1);
for k = 1:T
    % d = d_max * (2*rand(n,1) - 1); % New random disturbance
    d = 0;
    % subtract equilibrium from state to get nonzero equilibria
    u_cl = 0 + K*(x_cl(:,k) - [0 0 10 0]');
    u_cl = max(min(u_cl, 1), -1); % Saturate control input
    U1(:, k) = u_cl;
    x_cl(:,k+1) = A*x_cl(:,k) + B*u_cl + d;
end

% Plot results

subplot(4,1,1);
plot(1:T, X1);
title('System sampling');
xlabel('Time step');
ylabel('State');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(4,1,2);
plot(1:T, D0);
title('Disturbance signal');
xlabel('Time step');
ylabel('disturbance');
legend('d1', 'd2', 'd3', 'd4');
grid on;

subplot(4,1,3);
plot(0:T, x_cl);
title('Closed-loop System Response');
xlabel('Time step');
ylabel('State');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(4,1,4);
plot(1:T, U1);
title('Closed-loop Control Input');
xlabel('Time step');
ylabel('Input');
legend('u');
grid on;



% Analyze closed-loop eigenvalues
cl_eig = eig(A - B*K);
disp('Closed-loop eigenvalues:');
disp(cl_eig);