%{
furuta_verify.m
Holden Martinez
9/7/2026
Loads the derived model and checks energy/angular momentum conservation
across several initial conditions, with zero torque.
%}

load('furuta_model.mat', 'M_func', 'f_func');
p = furuta_params();
tau_fn = @(t,x) [0;0];   % zero torque throughout

tspan = [0 10];
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

x0_list = {
    [0; 0.3; 0;   0],    % rest start
    [0; 0.3; 1.5; 0],    % nonzero thetadot -- exercises Coriolis terms
    [0; 0.3; 1.0; 0.5]   % nonzero thetadot and alphadot
};
labels = {'Rest start', 'Nonzero thetadot', 'Nonzero thetadot and alphadot'};

figure;
for k = 1:length(x0_list)
    x0 = x0_list{k};

    [t_out, x_out] = ode45(@(t,x) furuta_ode(t, x, M_func, f_func, tau_fn, p), ...
        tspan, x0, opts);

    n = size(x_out, 1);
    E_vals = zeros(n,1);
    p_vals = zeros(n,1);

    for i = 1:n
        E_vals(i) = E(x_out(i,:), M_func, p);
        p_vals(i) = p_theta(x_out(i,:), M_func, p);
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