function V = j2_potential(r, p)
r = r(:);
rn = norm(r);
V = -p.mu/rn + p.mu*p.J2*p.Re^2*(3*r(3)^2/rn^2 - 1)/(2*rn^3);
end
