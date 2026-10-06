function x = mekf_init(qMeas, wMeas, p)
% Packed filter state x = [qHat(4); wHat(3); bHat(3); dHat(3); P(144)] started from the first
% star tracker fix: attitude from the measurement, rate from the first gyro sample, zero bias and
% disturbance torque estimates, and a diagonal covariance of the error state
% [attitude(3); rate(3); gyro bias(3); disturbance torque(3)].
P0 = blkdiag(p.starSigma^2*eye(3), (p.gyroArw^2/p.gyroStep)*eye(3), ...
             p.mekfBiasSigma0^2*eye(3), p.mekfTorqueSigma0^2*eye(3));
q0 = qMeas(:)/norm(qMeas);
x = [q0; wMeas(:); zeros(3, 1); zeros(3, 1); P0(:)];
end
