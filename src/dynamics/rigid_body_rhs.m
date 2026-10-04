function dx = rigid_body_rhs(x, inertia, torque)
% State x = [quaternion (scalar first, inertial to body); body rate], torque in the body frame.
q = x(1:4);
w = x(5:7);
dx = [quat_rate(q, w); inertia \ (torque(:) - cross(w, inertia*w))];
end
