%{
furuta_params.m
Holden Martinez
9/7/2026
Store parameters for easy future updating and reference, returns struct
%}
function p = furuta_params()
    p.L_a = 0.2; p.l_a = 0.1; p.l_p = 0.15;
    p.m_a = 0.05; p.m_p = 0.02; p.g = 9.81;
end