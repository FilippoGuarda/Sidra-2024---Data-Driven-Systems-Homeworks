
clear all
close all
rng(42)
%% system
n = 4; % dimension state, which we assume to be known
m = 1; % dimension input
p = 1; % dimension output
T = 100; % number of samples
J = T-n;
gam = 1;
A = [0.5780 0.8492 0.4220 0.1508;
     -0.6985 0.5780 0.6985 0.4220;
     0.4220 0.1508 0.5780 0.8492;
     0.6985 0.4220 -0.6985 0.5780];
B = [0.4610; 0.8492; 0.0390; 0.1508];
C = [0 0 1 0];
D = 0;  % Assuming D = 0 as in the original code

% Get system dimensions
[n, m] = size(B);  % n: number of states, m: number of inputs

%% Minimality test
OB = obsv(A,C);
if rank(OB) < n
    disp('System not observable');
    return
end

CO = ctrb(A,B);
if rank(CO) < n
    disp('System not controllable');
    return
end

%% I/O representation
[num,den] = ss2tf(A,B,C,D);
A_coeff = den(2:end);  % Remove leading 1
B_coeff = num(2:end);  % Remove leading 0
A_coeff = fliplr(A_coeff);
B_coeff = fliplr(B_coeff);

%% Auxiliary system
A_cal = zeros(2*n,2*n);
app2 = eye(n-1);
A_cal(1:n-1,2:n) = app2;
A_cal(n+1:2*n-1,n+2:2*n) = app2;
A_cal(n,:) = [-A_coeff B_coeff];

F_cal = zeros(2*n,2*n);
F_cal(1:n-1,2:n) = app2;
F_cal(n+1:2*n-1,n+2:2*n) = app2;

B_big = zeros(2*n,1);
B_big(end,1) = 1;
B_cal = B_big;

C_big = [-A_coeff B_coeff];

L_big = zeros(2*n,1);
L_big(n,1) = 1;
L_cal = L_big;

A_cal2 = F_cal + L_cal*[-A_coeff B_coeff];  % same as A_cal above

CO_big_sys = ctrb(A_cal,B_cal);
if rank(CO_big_sys) < 2*n
    disp('System A_cal, B_cal not reachable');
    return
end
%% data acquisition
X = zeros(n,T); % storage, corresponds to X_{0,T}
U = (2*gam).*rand(m,T+1)-gam; % storage, corresponds to U_{0,1,T}
Y = zeros(m,T); % storage, corresponds to Y_{0,1,T}
x = (2*gam).*rand(n,1)-gam; % initial conditions
for i =1:T+1
    u = U(:, i);
    X(:,i) = x;
    Y(:,i) = C*x;
    x = A*x+B*u;
end

M = zeros(2*n,J+1); % to construct matrices Phi0, Phi1
for i =1:n
    M(i,:) = Y(1,i:i+J);
end

for i =1:n
    M(n+i,:) = U(1,i:i+J);
end
Phi0 = M;
U0 = U(1,n+1:n+J+1);
N = [U0;Phi0];
if rank(N) < 2*n+1
    disp('PE condition failed');
    return
end

Phi_aux = [Y(1,J+2:J+n+1)'; U(1,J+2:J+n+1)'];
Phi1 = [Phi0(:,2:end) Phi_aux];
%% test on the identity A_cal*Phi0+B_cal*U0 = Phi1
if norm(A_cal*Phi0+B_cal*U0 - Phi1) > 1e-5
    disp('numerical problems');
    return
end




%% controller design (using CVX)
cvx_begin sdp
    variable Q(J+1,2*n)
    variable P(2*n,2*n) symmetric
    [P-eye(2*n) Phi1*Q; (Phi1*Q)' P] >= 0;
    P==Phi0*Q;
cvx_end

if strcmpi(cvx_status, 'Solved')
    disp('Robust controller design successful');
else
    error('Robust controller design failed');
end 

K_cal = U0*Q/P;
A_closed_loop_aux=A_cal+B_cal*K_cal;
disp('Aux system closed-loop eigenvalues (modulus)'); disp(abs(eig(A_closed_loop_aux)));
A_closed_loop=[A B*K_cal; L_cal*C F_cal+B_cal*K_cal];
disp('Closed-loop system eigenvalues (modulus)'); disp(abs(eig(A_closed_loop)));

% Simulate closed-loop system
x_cl = zeros(n, T+1);
U1 = zeros(1, T);
x_cl(:,1) = x(:,1); % Use the same initial condition
u_ctl = zeros(n, T+1);
Y = zeros(1, T);
Xi = zeros(2*n, 2*n);
Xe = [0 0 100 0]'; 
for k = 1:T
    x_cl(:, k+1) = A*x_cl(:, k) + B*K_cal*Xi(:, k);
    Xi(:, k+1) = L_cal*C*(x_cl(:, k) - Xe) + (F_cal+B_cal*K_cal)*Xi(:, k);
    Y(:, k) = C*x_cl(:, k);
    U1(:, k) = K_cal*Xi(:, k);
    U1(:, k) = max(min(U1(:, k), 1), -1);
end




% Plot results

subplot(4,1,1);
plot(0:T, X);
title('System sampling');
xlabel('Time step');
ylabel('State');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(4,1,2);
plot(0:T, x_cl);
title('Output feedback states');
xlabel('Time step');
ylabel('State');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(4,1,3);
plot(1:T, Y);
title('Closed-loop System Response');
xlabel('Time step');
ylabel('Output');
legend('Y');
grid on;

subplot(4,1,4);
plot(1:T, U1);
title('Closed-loop Control Input');
xlabel('Time step');
ylabel('Input');
legend('u');
grid on;