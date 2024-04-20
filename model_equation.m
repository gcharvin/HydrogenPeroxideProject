function dYdt = model_equation(t, Y, param)
    H = Y(1); A = Y(2); R = Y(3); Q = Y(4); S = Y(5);
    % Unpack parameters
    gamma = param.gamma; eta = param.eta; Hext = param.Hext_max;
    ep = param.ep; alpha = param.alpha; Ka = param.Ka;
    v = param.v; d = param.d; fa_max = param.fa_max; fr_max = param.fr_max;
    Htox = param.Htox; n = param.n; m = param.m;
    Ac = param.Ac; kp = param.kp; km = param.km; mu0 = param.mu0; mA = param.mA; mR = param.mR;

     fa = fa_max * H / (H + Ka);
    fr = fr_max;
    Kh = gamma / eta;
    Aoxi = H / (Kh + H) * A;
    mu = growth_rate(v, R, fr, mu0, mR, mA);
    G = kp * R * Aoxi^m / (Aoxi^m + Ac^m) - km * Q * Ac^m / (Aoxi^m + Ac^m);
    Jdamage = d * H^n / (H^n + Htox^n) * H * H * R;
    
    dHdt = ep + alpha * (Hext - H) - gamma * A * H / (Kh + H);
    dAdt = v * fa * R / mA - mu * A;
    dRdt = v * fr * R / mR - mu * R - Jdamage - G;
    dQdt = G - mu * Q;
    dSdt = Jdamage - mu * S;
    
    dYdt = [dHdt; dAdt; dRdt; dQdt; dSdt];
end

% % Example of running a simulation
% param = get_params();
% tspan = [0 100];
% Y0 = [0, 0.1, 1, 0.1, 0]; % Initial conditions
% [t, Y] = ode15s(@(t, Y) model_equation(t, Y, param), tspan, Y0);


