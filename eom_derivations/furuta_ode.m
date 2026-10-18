%{
furuta_ode.m
Holden Martinez
9/7/2026
%}
function dxdt = furuta_ode(t, x, M_func, f_func, tau_fn, p)
    theta = x(1); alpha = x(2); thetadot = x(3); alphadot = x(4);
    
    M_num = M_func(p.L_a, alpha, p.l_p, p.m_a, p.m_p, theta);
    f_num = f_func(p.L_a, alpha, alphadot, p.g, p.l_p, p.m_p, theta, thetadot);
    
    tau = tau_fn(t, x);     % <-- evaluated fresh every call, can depend on t and/or x
    qddot = M_num \ (tau + f_num);
    
    dxdt = [thetadot; alphadot; qddot(1); qddot(2)];
end