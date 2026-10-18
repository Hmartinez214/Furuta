%furuta_derive.m
% 10/17/2026
% Holden Martinez
% File intended to verify that algebra done by hand is equivalent for
% Furuta Pendulum (using Euler-Lagrange)
% Saves model in furuta_model.mat


syms t
syms L_a l_p h m_p m_a g real
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

%v_psquared_hand= L_a^2*thetadot^2-2*L_a*l_p*alphadot*thetadot*cos(alpha)+thetadot^2*l_p^2*sin(alpha)^2+l_p^2*alphadot^2;

%Check the difference between matlab and hand version
% 1 -> equivalent 
% 0 -> invalid

%diff_expr = simplify(v_psquaredclean - v_psquared_hand, 'Steps', 100);
%simplify(diff_expr, 'IgnoreAnalyticConstraints', true)
%isAlways(diff_expr == 0)

%Moments of Inertia
I_a= (1/3)*m_a*L_a^2;
I_p= (1/3)*m_p*l_p^2;

%Kinetic
k_a=(1/2)*I_a*Dtheta^2;
k_prot=(1/2)*I_p*(Dalpha^2+Dtheta^2*sin(alpha)^2);
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
%thetamotion_hand = (1/3)*m_a*L_a^2*thetadoubledot+m_p*L_a^2*thetadoubledot-m_p*L_a*l_p*alphadoubledot*cos(alpha)+m_p*L_a*l_p*alphadot^2*sin(alpha)+m_p*thetadoubledot*l_p^2*sin(alpha)^2+m_p*thetadot*2*sin(alpha)*cos(alpha)*alphadot*l_p^2;
%check if difference
%diff_theta_eom = simplify(thetamotionclean - thetamotion_hand, 'Steps', 200);
%isAlways(diff_theta_eom == 0)

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
%alphamotion_hand=(1/12)*m_p*l_p^2*alphadoubledot-m_p*l_p*L_a*thetadoubledot*cos(alpha)+m_p*L_a*l_p*thetadot*alphadot*sin(alpha)+m_p*l_p^2*alphadoubledot-m_p*L_a*l_p*alphadot*thetadot*sin(alpha)-m_p*thetadot^2*l_p^2*cos(alpha)*sin(alpha)-m_p*g*l_p*sin(alpha);
%check if different
%diff_alpha_eom = simplify(alphamotionclean - alphamotion_hand, 'Steps', 200);
%isAlways(diff_alpha_eom == 0)

thetamotion_numready = subs(thetamotionclean, [theta, alpha], [theta_s, alpha_s]);
alphamotion_numready = subs(alphamotionclean, [theta, alpha], [theta_s, alpha_s]);


eqns= [thetamotion_numready == 0, alphamotion_numready ==0]; %EL equations
vars= [thetadoubledot, alphadoubledot];
[M,f]= equationsToMatrix(eqns, vars); %combine into matrix
%convert M, F, into numeric functions
M_func = matlabFunction(M, 'Vars', [L_a, alpha_s, l_p, m_a, m_p, theta_s]);
f_func = matlabFunction(f, 'Vars', [L_a, alpha_s, alphadot, g, l_p, m_p, theta_s, thetadot]);

save('furuta_model.mat', 'M_func', 'f_func');
save('furuta_symbolic_model.mat', 'thetamotion_numready', 'alphamotion_numready');