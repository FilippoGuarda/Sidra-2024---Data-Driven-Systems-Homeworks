%% NUMERICAL EXAMPLE ON CONTRACTIVITY
%  system: one-link robot arm
%  REMARKS
%  In this example we pick Q(x)=cox(x1) or Q(x)=[cox(x1) x1^2 sin(x2)]'

%% Initialization
clear,clc
rng(1);

%% System parameters
global K
global A
global B

n = 4; % system dimensions
m = 1;
T = 12; % number of samples
% system parameters
A = [0 1 0 0; ... % precomputed A matrix
    -2 -0.75 1 0; ...
    0 0 0 1; ...
    -1.3 0 0.66 -0.2];
B = [0; 0; 0; 6.6]; % precomputed B matrix
mag = 0.1;
rank(ctrb(A, B))

T    = 10;  % number of samples
Ts   = 0.1; % sampling interval
Tsim = T*Ts; % duration of simulation

mag = 0.1; % magnitude of initial conditions
x0  = (2*mag).*rand(n,1)-mag; % initial state

sim('data_collection_arm');

x  = state.signals.values'; 
xd = state_deriv.signals.values'; 
u  = input.signals.values'; 

X0  = x(:,1:T);
U0  = u(:,1:T);
X1  = xd(:,1:T);

s = 5; % dimension of Z(x)
Z0  = [X0;cos(X0(1,:))];
RQ = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % Bound for Jacobian of Q(x)

rank([U0; Z0]) % full row rank check

r = n; % column number of RQ

cvx_begin sdp
    variable P1(n,n) symmetric
    variable Y1(T,n)
    variable G2(T,s-n)
    variable a 
    P1 >= 0*eye(n);
    a >= 0;
    Z0*Y1 == [P1;zeros(s-n,n)];
    [X1*Y1+transpose(X1*Y1)+a*eye(n) X1*G2 P1*RQ;
        transpose(X1*G2) -eye(s-n) zeros(s-n,r);
        transpose(P1*RQ) zeros(r,s-n) -eye(r)] <= 0;
    Z0*G2 == [zeros(n,s-n);eye(s-n)];
cvx_end

G1 = Y1/P1;
G  = [G1 G2];
K  = U0*G; 
closed = X1*G;

mx0 = 1;   % magnitude of initial conditions
x0  = (2*mx0).*rand(n,1)-mx0;
tspan = [0,100];  % duration of simulation

[t,x] = ode45(@arm,tspan,x0);

x_star = x(end,:)';
u_star = K*[x_star;cos(x_star(1))]; % equlibrium point of closed-loop system

figure
plot(t,x(:,1),'r','LineWidth',1);
hold on;
plot(t,x(:,2),'b','LineWidth',1);
hold on;
plot(t,x(:,3),'g','LineWidth',1);
hold on;
plot(t,x(:,4),'k','LineWidth',1);
xlabel('t');
legend('x(1) ','x(2)','x(3) ','x(4)');

function dxdt = arm(t,x)  
    global K
    global A
    global B

    dxdt = zeros(4,1);
    u = K*[x;cos(x(1))];
    dxdt = A*x + B*u;
    dxdt(2) = dxdt(2) - 1.96*cos(x(1));
    
end

% % Simulate closed-loop system
% x_cl = zeros(n, T+1);
% U1 = zeros(1, T);
% x_cl(:,1) = x(:,1); % Use the same initial condition
% u_ctl = zeros(n, T+1);
% Y = zeros(1, T);
% Xi = zeros(2*n, 2*n);
% for k = 1:T
%     x_cl(:, k+1) = A*x_cl(:, k) + B*K_cal*Xi(:, k);
%     Xi(:, k+1) = L_cal*C*x_cl(:, k) + (F_cal+B_cal*K_cal)*Xi(:, k);
%     Y(:, k) = C*x_cl(:, k);
%     U1(:, k) = K_cal*Xi(:, k);
% end
