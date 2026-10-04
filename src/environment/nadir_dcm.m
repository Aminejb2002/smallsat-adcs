function A = nadir_dcm(r, v)
% Inertial-to-body DCM for nadir pointing: z nadir, y negative orbit normal, x along track.
r = r(:);
v = v(:);
z = -r/norm(r);
h = cross(r, v);
y = -h/norm(h);
x = cross(y, z);
A = [x'; y'; z'];
end
