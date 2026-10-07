function [tau, parts] = disturbance_block(t, r, v, q, p)
% Disturbance torques in the body frame: gravity gradient, solar pressure, drag, residual magnetic dipole.
% parts holds one column per source in that order, tau is their sum.
A = quat_to_dcm(q);
sun = sun_direction(p.epochDays + t/86400);
lit = ~in_shadow(r, sun, p.Re);
parts = [torque_gravity_gradient(A*r(:), p.inertiaTruth, norm(r), p), ...
         torque_srp(A, sun, lit, p), ...
         torque_drag(A, r, v, p), ...
         cross(p.residualDipole, A*magnetic_field_eci(r, t, p))];
tau = sum(parts, 2);
end
