function u = pointing_law(q, w, qref, wref, p)
% PD law on the attitude error quaternion with a slew-rate limit. u is the torque wanted on the
% body; the wheels must apply -u. q is inertial to body, wref is given in the reference frame.
s = size_pointing(p);
qe = quat_mul([qref(1); -qref(2:4)], q(:));
if qe(1) < 0
    sgn = -1;
else
    sgn = 1;
end
werr = w(:) - quat_to_dcm(qe)*wref(:);
theta = 2*sgn*qe(2:4);
wdes = -min(max(s.wn*theta/(2*s.zeta), -s.wSlew), s.wSlew);
u = -s.kd.*(werr - wdes);
end
