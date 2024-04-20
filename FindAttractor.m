function [t, Y, converged] = FindAttractor(param)
    Tol = 1e-4;
    tmax = 1e2;
    while tmax < 1e8
        tlist = logspace(-2, log10(tmax / 10.0), 64);
        tlist = [tlist, tmax * 0.99];
        opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
        [t, Y] = ode15s(@(t, Y) model_equation(t, Y, param), tlist, [0, 0.1, 1, 0.1, 0], opts);
        x = Y(NearestTimePoint(t, tmax / 10), :);
        y = Y(end, :);
        
        dist = max(abs(x - y) ./ (x + y));
        if dist < Tol
            converged = true;
            break;
        else
            tmax = tmax * 2.0;
            converged = false;
        end
    end
end


function index = NearestTimePoint(tlist, t)
    [~, index] = min(abs(tlist - t));
end
