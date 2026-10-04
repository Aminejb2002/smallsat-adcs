function s = sun_direction(daysJ2000)
% Low-precision Sun direction in the inertial frame, days since J2000.
d2r = pi/180;
L = (280.460 + 0.9856474*daysJ2000)*d2r;
g = (357.528 + 0.9856003*daysJ2000)*d2r;
lam = L + 1.915*d2r*sin(g) + 0.020*d2r*sin(2*g);
obl = (23.4393 - 3.563e-7*daysJ2000)*d2r;
s = [cos(lam); sin(lam)*cos(obl); sin(lam)*sin(obl)];
end
