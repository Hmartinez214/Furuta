%{
furuta_lqr.m
Holden Martinez
9/10/2026
Jacobian and Linearization process for the eoms of furuta pendulum
%}  
clear
close all
load('furuta_model.mat', 'M_func', 'f_func');
load('furuta_symbolic_model.mat','alphamotion_numready','thetamotion_numready');
p = furuta_params();
syms theta_s alpha_s real
syms thetadot alphadot thetadoubledot alphadoubledot
syms tau real
syms L_a l_a l_p m_a m_p g real

%Tuning parameters
R_o=10e3;
alpha_o=.633;
t_o=5;

B = [1; 0];
q = [theta_s; alpha_s];
qdot = [thetadot; alphadot];

eqns = [thetamotion_numready == B(1)*tau, alphamotion_numready == B(2)*tau];
vars = [thetadoubledot, alphadoubledot];
[M_sym, f_sym] = equationsToMatrix(eqns, vars);   % note: now WITH tau in eqns

qddot_sym = M_sym \ f_sym;   % symbolic M^{-1} * f

x_sym = [theta_s; alpha_s; thetadot; alphadot];
f_full = [thetadot; alphadot; qddot_sym(1); qddot_sym(2)];

A=jacobian(f_full,x_sym);
x_eq = [0, 0, 0, 0];   % [theta, alpha, thetadot, alphadot]
tau_eq = 0;
A_eq = subs(A, [theta_s, alpha_s, thetadot, alphadot, tau, L_a, l_a, l_p, m_a, m_p, g], ...
               [x_eq, tau_eq, p.L_a, p.l_a, p.l_p, p.m_a, p.m_p, p.g]);

A_eq = double(A_eq);   % convert from symbolic to plain numeric matrix
Bvec = jacobian(f_full, tau);
B_eq = double(subs(Bvec, [theta_s, alpha_s, thetadot, alphadot, tau, L_a, l_a, l_p, m_a, m_p, g], ...
                            [x_eq, tau_eq, p.L_a, p.l_a, p.l_p, p.m_a, p.m_p, p.g]));
eigvals=eig(A_eq);

rank(ctrb(A_eq, B_eq))
Q = diag([1 10 0.1 0.1]); %theta alpha and respective velocities
R = R_o;
K = lqr(A_eq, B_eq, Q, R);
   tau_fn = @(t,x) [-K*x; 0];
   [t_out, x_out] = ode45(@(t,x) furuta_ode(t,x,M_func,f_func,tau_fn,p), [0 t_o], [0; alpha_o; 0; 0]);

plot(t_out, x_out(:,1:2));
legend('\theta','\alpha');
xlabel('time(s)');
ylabel('Angles (rad)')
figure;
plot(t_out, x_out(:,3:4))
legend('$\dot{\theta}$','$\dot{\alpha}$', 'Interpreter', 'latex');
xlabel('time(s)');
ylabel('Angular Velocity (rad/s)')
Poles = eig(A_eq - B_eq*K)
u_matrix= -K.*x_out; %compute and plot torque
u=sum(u_matrix,2); 
figure;
plot(t_out, u);
legend('torque');
xlabel('time(s)');
ylabel('torque(Nm)')


