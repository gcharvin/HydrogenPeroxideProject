function sol = ComputeODE(param, steady, tlist)
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    sol = ode15s(@(t, Y) model_equation(t, Y, param), [0, max(tlist)], steady(end, :), opts);
end
