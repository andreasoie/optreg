% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Fl�ten
addpath(genpath("utils"))

%% Initialization and model definition
init08; % Change this to the init file corresponding to your helicopter

% Constants
DELTA_T	= 0.25; % sampling time
ALPHA = 0.2;
BETA = 20;
LAMBDA_T = (2*pi) / 3;

% Parameters
q1 = 1;
q2 = 1;

% Discrete time system model. x = [lambda r p p_dot e e_dot]'
A1 = [1 DELTA_T 0 0 0 0;
      0 1 -K_2*DELTA_T 0 0 0;
      0 0 1 DELTA_T 0 0;
      0 0 -K_1*K_pp*DELTA_T -K_1*K_pd*DELTA_T+1 0 0;
      0 0 0 0 1 DELTA_T;
      0 0 0 0 -K_3*K_ep*DELTA_T -K_3*K_ed*DELTA_T+1];
B1 = [0 0;
      0 0;
      0 0;
      K_1*K_pp*DELTA_T 0;
      0 0;
      0 K_3*K_ep*DELTA_T];

% Number of states and inputs
mx = size(A1, 2); % Number of states (number of columns in A)
mu = size(B1, 2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % p_dot
x6_0 = 0;                               % p_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';            % Initial values

% Time horizon and initialization
N  = 40;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

MAX_ANGLE = 30;
MAX_RADIAN = (MAX_ANGLE*pi) / 180;

% Bounds
ul 	    = -MAX_RADIAN;                   % Lower bound on control
uu 	    =  MAX_RADIAN;                    % Upper bound on control

uLower = [ul; -inf];
uUpper = [uu;  inf];

xl(1:mx, 1) = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu(1:mx, 1) = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb, vub]      = gen_constraints(N, M, xl, xu, uLower, uUpper);
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx, mx);
Q1(1,1) = 2;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;                            % Weight on state x5
Q1(6,6) = 0;                            % Weight on state x6

P1 = zeros(2, 2);                      
P1(1, 1) = q1;                          % Weight on input "u1" pitch
P1(2, 2) = q2;                          % Weight on input "u2" elevation

Q = 2*gen_q(Q1, P1, N, M);   % Generate Q, hint: gen_q
c = zeros(N*mx+M*mu, 1);     % Generate c, this is the linear constant term in the QP

%% Generate system matrixes for linear model
Aeq = gen_aeq(A1, B1, N, mx, mu);             % Generate A, hint: gen_aeq
beq = zeros(size(Aeq, 1), 1);
beq(1:mx) = A1*x0;

% Inequalities
Aineq = [];
bineq = [];

% Minimize
CF = @(z) (1/2) * (z') * Q * z;

opt = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunEvals', 10e4);

tic
[Z, ZVAL, EXITFLAG] = fmincon(CF, z0, Aineq, bineq, Aeq, beq, vlb, vub, @c_constrain, opt);
toc

% Extract signals
x_values = Z(1:N*mx, 1);
u_values = Z(N*mx+1:end, 1);

% Allocate states
x1 = x_values(1:mx:end);  % State x1 from solution
x2 = x_values(2:mx:end);  % State x2 from solution
x3 = x_values(3:mx:end);  % State x3 from solution
x4 = x_values(4:mx:end);  % State x4 from solution
x5 = x_values(5:mx:end);  % State x5 from solution
x6 = x_values(6:mx:end);  % State x6 from solution

% Allocate inputs
u1  = u_values(1:mu:end); % Control input "u1" from solution
u2  = u_values(2:mu:end); % Control input "u2" from solution

num_variables = 5/DELTA_T;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];

x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

%% Plotting
t1 = 0:DELTA_T:DELTA_T*(length(u1)-1);
t2 = 0:DELTA_T:DELTA_T*(length(u2)-1);
assert(size(t1, 2) == size(t1, 2), 'DimentionException')

% figure(1337)
% 
% subplot(811)
% stairs(t1, u1, "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Input: Pitch")
% 
% subplot(812)
% stairs(t2, u2, "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Input: Elevation")
% 
% subplot(813)
% plot(t1, x1,"m", t1, x1, "*", "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Travel")
% 
% subplot(814)
% plot(t1, x2, "m", t1, transpose(x2), "*", "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Travel Rate")
% 
% subplot(815)
% plot(t1, x3, "m", t1, transpose(x3), "*", "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Pitch")
% 
% subplot(816)
% plot(t1, x4,"m", t1, transpose(x4),"*", "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Pitch Rate")
% 
% subplot(817)
% plot(t1, x5, "m", t1, transpose(x5),"*", "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Elevation")
% 
% subplot(818)
% plot(t1, x6, "m", t1, transpose(x6),"*", "LineWidth", 2)
% grid
% xlabel("time (s)")
% ylabel("radians (rad)")
% title("Elevation Rate")

%% LQ Controller

q_travel = 150;
q_travel_dot = 1;
q_pitch = 0.1;
q_pitch_dot = 15;
q_elevation = 10;
q_elevation_dot = 20;

r_pitch = 0.1;
r_elevation = 0.1;

Q1 = diag([q_travel q_travel_dot q_pitch q_pitch_dot q_elevation q_elevation_dot]); 
R1 = diag([r_pitch r_elevation]);

% Calculates the optimal gain matrix K for state-space
[K_opt, S, CLP] = dlqr(A1, B1, Q1, R1);

% Save stuff for Simulink aswell as plotting options

% Pitch
Up_opt= timeseries(u1, t1); 
% Elevation
Ue_opt= timeseries(u2, t2); 

X_opt = timeseries([x1 x2 x3 x4 x5 x6], t1);

T_opt = t1;
