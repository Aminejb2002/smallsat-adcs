function x = mekf_init(qMeas, p)
% Packed filter state x = [qHat(4); bHat(3); P(36)] started from the first star tracker fix,
% zero bias estimate and a diagonal covariance.
P0 = blkdiag(p.starSigma^2*eye(3), p.mekfBiasSigma0^2*eye(3));
q0 = qMeas(:)/norm(qMeas);
x = [q0; zeros(3, 1); P0(:)];
end
