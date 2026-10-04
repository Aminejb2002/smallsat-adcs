function ic = initial_state(p, pitchOffsetDeg)
% Nadir-pointing start on the dawn-dusk Sun-synchronous orbit, body rate equal to the orbit rate.
% Fields r, v, q, w are column vectors; q is scalar first, inertial to body.
if nargin < 2
    pitchOffsetDeg = 2;
end
a = p.Re + p.altitude;
n = sqrt(p.mu/a^3);
inc = sso_inclination(a, p);
raan = dawn_dusk_raan(p.epochDays);
[r0, v0] = circular_state(a, inc, raan, 0, p);
th = pitchOffsetDeg*pi/180;
ry = [cos(th) 0 -sin(th); 0 1 0; sin(th) 0 cos(th)];
A0 = ry*nadir_dcm(r0, v0);
hHat = cross(r0, v0)/norm(cross(r0, v0));
ic.r = r0;
ic.v = v0;
ic.q = dcm_to_quat(A0);
ic.w = A0*(n*hHat);
end
