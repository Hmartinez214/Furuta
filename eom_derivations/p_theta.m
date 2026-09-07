%{
p_theta.m
Holden Martinez
9/7/2026
Computes angular momentum conjugate to theta
Theta is cyclic, p_theta should be constant
%}
function P_theta = p_theta(x, M_func, p)
    P_theta = [1 0] * M_func(p.L_a, x(2), p.l_a, p.l_p, p.m_a, p.m_p, 0) * [x(3); x(4)];
end