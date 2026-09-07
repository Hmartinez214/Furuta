%Pendulum_verification_correct.m
% 9/7/2026
% Holden Martinez
% File intended to verify that algebra done by hand is equivalent for
% Furuta Pendulum (using Euler-Lagrange)


syms t
syms L_a l_a L_p l_p h m_p m_a g real
syms theta(t) alpha(t)

syms alphadot alphadoubledot
syms thetadot thetadoubledot

syms theta_s alpha_s real

Dtheta = diff(theta, t);   
Dalpha = diff(alpha, t);

%position

x_p= L_a*sin(theta)-l_p*sin(alpha)*cos(theta);
y_p= L_a*cos(theta)+l_p*sin(alpha)*sin(theta);
z_p= h+l_p*cos(alpha);

%velocity 

x_pdot = diff(x_p, t);
y_pdot= diff(y_p, t);
z_pdot= diff(z_p, t);
%combined and squared
v_psquared=x_pdot^2+y_pdot^2+z_pdot^2;
v_psquaredclean=subs(v_psquared, [diff(alpha, t), diff(theta, t)], [alphadot, thetadot]);

v_psquared_hand= L_a^2*thetadot^2-2*L_a*l_p*alphadot*thetadot*cos(alpha)+thetadot^2*l_p^2*sin(alpha)^2+l_p^2*alphadot^2;

%Check the difference between matlab and hand version
% 1 -> equivalent 
% 0 -> invalid

diff_expr = simplify(v_psquaredclean - v_psquared_hand, 'Steps', 100);
%simplify(diff_expr, 'IgnoreAnalyticConstraints', true)
isAlways(diff_expr == 0)

%Moments of Inertia
I_a= (1/3)*m_a*l_a^2;
I_p= (1/12)*m_p*l_p^2;

%Kinetic
k_a=(1/2)*I_a*Dtheta^2;
k_prot=(1/2)*I_p*Dalpha^2;
k_ptrans=(1/2)*m_p*v_psquared;
T=k_a+k_prot+k_ptrans;
%Potential
V=m_p*g*l_p*cos(alpha);

L=T-V;
%EOM -- wrt theta
%calculate
thetamotion=diff(diff(L, Dtheta),t)-diff(L, theta);
%substitute to make look cleaner
thetamotionclean = subs(thetamotion, ...
    [diff(theta,t,2), diff(alpha,t,2), diff(theta,t), diff(alpha,t)], ...
    [thetadoubledot, alphadoubledot, thetadot, alphadot]); 
%simplify
thetamotionclean_simplified = simplify(thetamotionclean, 'Steps', 200);
%enter value calculated by hand
thetamotion_hand = (1/3)*m_a*l_a^2*thetadoubledot+m_p*L_a^2*thetadoubledot-m_p*L_a*l_p*alphadoubledot*cos(alpha)+m_p*L_a*l_p*alphadot^2*sin(alpha)+m_p*thetadoubledot*l_p^2*sin(alpha)^2+m_p*thetadot*2*sin(alpha)*cos(alpha)*alphadot*l_p^2;
%check if difference
diff_theta_eom = simplify(thetamotionclean - thetamotion_hand, 'Steps', 200);
isAlways(diff_theta_eom == 0)

%EOM -- wrt alpha
%calculate
alphamotion=diff(diff(L, Dalpha),t)-diff(L,alpha);
%subs to make cleaner
alphamotionclean=subs(alphamotion, ...
    [diff(theta,t,2), diff(alpha,t,2), diff(theta,t), diff(alpha,t)], ...
    [thetadoubledot, alphadoubledot, thetadot, alphadot]); 
%simplify
alphamotionclean_simplified = simplify(alphamotionclean, 'Steps', 200);
%enter value calculated by hand
alphamotion_hand=(1/12)*m_p*l_p^2*alphadoubledot-m_p*l_p*L_a*thetadoubledot*cos(alpha)+m_p*L_a*l_p*thetadot*alphadot*sin(alpha)+m_p*l_p^2*alphadoubledot-m_p*L_a*l_p*alphadot*thetadot*sin(alpha)-m_p*thetadot^2*l_p^2*cos(alpha)*sin(alpha)-m_p*g*l_p*sin(alpha);
%check if different
diff_alpha_eom = simplify(alphamotionclean - alphamotion_hand, 'Steps', 200);
isAlways(diff_alpha_eom == 0)

thetamotion_numready = subs(thetamotionclean, [theta, alpha], [theta_s, alpha_s]);
alphamotion_numready = subs(alphamotionclean, [theta, alpha], [theta_s, alpha_s]);


eqns= [thetamotion_numready == 0, alphamotion_numready ==0]; %EL equations
vars= [thetadoubledot, alphadoubledot];
[M,f]= equationsToMatrix(eqns, vars); %combine into matrix
%convert M, F, into numeric functions
M_func = matlabFunction(M, 'Vars', [L_a, alpha_s, l_a, l_p, m_a, m_p, theta_s]);
f_func = matlabFunction(f, 'Vars', [L_a, alpha_s, alphadot, g, l_p, m_p, theta_s, thetadot]);

function dxdt = furuta_ode(t, x, M_func, f_func, tau, L_a, l_a, l_p, m_a, m_p, g)
    theta = x(1);
    alpha = x(2);
    thetadot = x(3);
    alphadot = x(4);
    
    M_num = M_func(L_a, alpha, l_a, l_p, m_a, m_p, theta);
    f_num = f_func(L_a, alpha, alphadot, g, l_p, m_p, theta, thetadot);
    
    qddot = M_num \ (tau + f_num);
    
    dxdt = [thetadot; alphadot; qddot(1); qddot(2)];
end


% Energy and angular momentum verification

tau_scalar = 0;   % zero torque throughout
B = [1; 0];
tspan = [0 10];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

% numeric parameter values -- placeholders currently
L_a_val = 0.2; l_a_val = 0.1; l_p_val = 0.15; m_a_val = 0.05; m_p_val = 0.02; g_val = 9.81;

% Energy and angular momentum as functions of state x = [theta; alpha; thetadot; alphadot]
E = @(x) 0.5*[x(3);x(4)]' * M_func(L_a_val,x(2),l_a_val,l_p_val,m_a_val,m_p_val,0) * [x(3);x(4)] ...
        + m_p_val*g_val*l_p_val*cos(x(2));

p_theta = @(x) [1 0] * M_func(L_a_val, x(2), l_a_val, l_p_val, m_a_val, m_p_val, 0) * [x(3); x(4)];

% Test cases: rest start, nonzero thetadot (exercises Coriolis terms), and both nonzero
x0_list = {
    [0; 0.3; 0;   0],    % rest start
    [0; 0.3; 1.5; 0],    % nonzero thetadot
    [0; 0.3; 1.0; 0.5]   % nonzero thetadot and alphadot
};
labels = {'Rest start', 'Nonzero thetadot', 'Nonzero thetadot and alphadot'};

figure;
for k = 1:length(x0_list)
    x0 = x0_list{k};

    [t_out, x_out] = ode45(@(t,x) furuta_ode(t, x, M_func, f_func, B*tau_scalar, ...
        L_a_val, l_a_val, l_p_val, m_a_val, m_p_val, g_val), tspan, x0, opts);

    n = size(x_out, 1);
    E_vals = zeros(n,1);
    p_vals = zeros(n,1);

    for i = 1:n
        E_vals(i) = E(x_out(i,:));
        p_vals(i) = p_theta(x_out(i,:));
    end

    E_drift = (E_vals - E_vals(1)) / abs(E_vals(1));
    p_drift = (p_vals - p_vals(1));

    subplot(2, length(x0_list), k);
    plot(t_out, E_drift);
    ylabel('Relative energy drift');
    title(['Energy: ', labels{k}]);

    subplot(2, length(x0_list), k + length(x0_list));
    plot(t_out, p_drift);
    ylabel('p_\theta drift');
    title(['Momentum: ', labels{k}]);
end