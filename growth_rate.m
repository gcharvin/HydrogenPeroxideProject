
function mu = growth_rate(v, R, fr, mu0, mR, mA)
    mu = mu0 * v * R * (fr / mR + (1 - fr) / mA);
end
