function rho = atmos_density(h, p)
% Exponential density model about the reference altitude.
rho = p.rhoRef*exp(-(h - p.altitude)/p.scaleHeight);
end
