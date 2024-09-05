clear all
rng(1);

% System parameters
global K

n = 4; % dimension of state
m = 1; % dimension of input
p = 1; % dimension of output
d_max = 0.1;

A = [0 1 0 0; ... % precomputed A matrix
    -2 -0.75 1 0; ...
    0 0 0 1; ...
    -1.3 0 0.66 -0.2];
B = [0; 0; 0; 6.6]; % precomputed B matrix

target = pi/3; % target reference

T    = 10;  % number of samples
Ts   = 0.1; % sampling interval
Tsim = T*Ts; % duration of simulation

mag = 0.1; % magnitude of initial conditions
x0  = (2*mag).*rand(n+p,1)-mag; % initial state

sim('robot_arm_pi2');

x  = state.signals.values'; 
xd = state_deriv.signals.values'; 
u  = input.signals.values'; 

X0  = x(:,1:T);
U0  = u(:,1:T);
Z1  = xd(:,1:T);

s = 5; % dimension of Z(x)
Z0  = [X0;cos(X0(1,:))];
RQ = [1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % Bound for Jacobian of Q(x)

M = ones(1,T);
rank([U0; Z0; M]) % full row rank check

r = n; % column number of RQ

cvx_begin sdp
    variable P1(n+p,n+p) symmetric
    variable Y1(T,n+p)
    variable G2(T,s-n)
    variable a
    P1 >= 0*eye(n+p);
    a >= 0;
    Z0*Y1 == [P1;zeros(s-n,n+p)];
    [Z1*Y1+transpose(Z1*Y1)+a*eye(n+p) Z1*G2 P1*RQ;
        transpose(Z1*G2) -eye(s-n) zeros(s-n,r);
        transpose(P1*RQ) zeros(r,s-n) -eye(r)] <= 0;
    Z0*G2 == [zeros(n+p,s-n);eye(s-n)];
    M*[Y1 G2] == zeros(1,s+p);
cvx_end

G  = [Y1/P1 G2];
K  = U0*G; 
closed = Z1*G;

mag = 1;
x0  = (2*mag).*rand(5,1)-mag;
tspan = [0,100]; 
r = pi/3;

[t,x_cl] = ode45(@arm_integral,tspan,[x0; r]);

subplot(2,1,1);
plot(0:T, x(1:4, :));
title('System sampling');
xlabel('Time step');
ylabel('State');
legend('x1', 'x2', 'x3', 'x4');
grid on;

subplot(2,1,2);
plot(t, x_cl(:, 1:4));
title('PI states');
xlabel('Time step');
ylabel('State');
legend('x1', 'x2', 'x3', 'x4');
grid on;

function dxdt = arm_integral(t,x)  
    global K
    A = [0 1 0 0; ... % precomputed A matrix
    -2 -0.75 1 0; ...
    0 0 0 1; ...
    -1.3 0 0.66 -0.2];
    B = [0; 0; 0; 6.6]; % precomputed B matrix
    target = pi/3;
    noise = 0.1 * (2*rand(4,1) - 1);

    dxdt = zeros(6,1);
    u = K*[x(1:5);cos(x(1))];
    dxdt(1:4) = A*x(1:4) + B*u + noise;
    dxdt(2) = dxdt(2) - 1.96*cos(x(1));
    dxdt(5) = x(1) - target;
    dxdt(6) = 0;
end