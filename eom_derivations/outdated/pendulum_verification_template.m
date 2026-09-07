%pendulum_Verification_template.m
% 9/6/2026
% Holden Martinez
% File intended to verify that algebra done by hand is equivalent for
% Furuta Pendulum (using Euler-Lagrange)


syms t
syms L_a l_a L_p l_p h m_p m_a g real
syms theta(t) alpha(t)

syms alphadot alphadoubledot
syms thetadot thetadoubledot

Dtheta = diff(theta, t);   
Dalpha = diff(alpha, t);

%position

x_p= 
y_p= 
z_p= 

%velocity 

x_pdot = diff(x_p, t);
y_pdot= diff(y_p, t);
z_pdot= diff(z_p, t);

%combined and squared

v_psquared=x_pdot.^2+y_pdot.^2+z_pdot.^2;

%substitute "diff" terms for alphadot and thetadot

v_psquaredclean=subs(v_psquared, [diff(alpha, t), diff(theta, t)], [alphadot, thetadot]);

%below is the value I got calculating by hand

v_psquaredcalc= (L_a + l_p*sin(alpha))^2 * thetadot^2 + l_p^2*alphadot^2;

%Check the difference between matlab and my version
% 1 -> equivalent 
% 0 -> invalid

diff_expr = simplify(v_psquaredclean - v_psquaredcalc, 'Steps', 100);
simplify(diff_expr, 'IgnoreAnalyticConstraints', true)
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
thetamotion_hand = 
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
alphamotion_hand=
%check if different
diff_alpha_eom = simplify(alphamotionclean - alphamotion_hand, 'Steps', 200);
isAlways(diff_alpha_eom == 0)