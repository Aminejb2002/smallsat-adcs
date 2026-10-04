function B = magnetic_field_eci(r, t, p)
% Tilted-dipole field in the inertial frame; the dipole turns with the Earth.
r = r(:);
rn = norm(r);
rh = r/rn;
lat = -80.7*pi/180;
lon = 107.3*pi/180;
mEcef = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)];
th = p.wEarth*t;
m = [cos(th)*mEcef(1) - sin(th)*mEcef(2); sin(th)*mEcef(1) + cos(th)*mEcef(2); mEcef(3)];
B = p.dipole/rn^3*(3*(m'*rh)*rh - m);
end
