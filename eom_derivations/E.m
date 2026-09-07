%{
E.m
Holden Martinez
9/7/2026
Function to calculate the energy of the pendulum
%}
function Energy = E(x, M_func, p)
    Energy = 0.5*[x(3);x(4)]' * M_func(p.L_a, x(2), p.l_a, p.l_p, p.m_a, p.m_p, 0) * [x(3);x(4)] ...
            + p.m_p*p.g*p.l_p*cos(x(2));
end