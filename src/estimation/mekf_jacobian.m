function F = mekf_jacobian(w, h, inertia)
% Continuous error-state matrix of the filter for [attitude(3); rate(3); bias(3); torque(3)]
% about the rate w and wheel momentum h. Attitude error is in the body frame; the rate follows
% inertia*wdot = d - tauW - w x (inertia*w + h) with d the disturbance torque.
w = w(:);
v = inertia*w + h(:);
F = zeros(12);
F(1:3, 1:3) = -skew_matrix(w);
F(1:3, 4:6) = eye(3);
F(4:6, 4:6) = inertia\(skew_matrix(v) - skew_matrix(w)*inertia);
F(4:6, 10:12) = inertia\eye(3);
end
