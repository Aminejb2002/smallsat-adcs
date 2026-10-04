function inc = sso_inclination(a, p)
% Inclination of a circular Sun-synchronous orbit with semi-major axis a.
n = sqrt(p.mu/a^3);
rate = 2*pi/p.yearSeconds;
inc = acos(-rate/(1.5*n*p.J2*(p.Re/a)^2));
end
