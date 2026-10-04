function a = j2_accel(r, p)
r = r(:);
rn = norm(r);
zz = (r(3)/rn)^2;
k = 1.5*p.J2*p.mu*p.Re^2/rn^5;
a = -p.mu*r/rn^3 + k*[r(1)*(5*zz - 1); r(2)*(5*zz - 1); r(3)*(5*zz - 3)];
end
