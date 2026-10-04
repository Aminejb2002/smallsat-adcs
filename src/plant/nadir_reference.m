function [qref, wref] = nadir_reference(r, v)
% Nadir-pointing reference: attitude quaternion (inertial to body) and its rate in the reference frame.
r = r(:);
v = v(:);
h = cross(r, v);
qref = dcm_to_quat(nadir_dcm(r, v));
wref = [0; -norm(h)/(r'*r); 0];
end
