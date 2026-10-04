function A = quat_to_dcm(q)
% Scalar-first quaternion to the inertial-to-body direction cosine matrix.
q = q(:);
qv = q(2:4);
A = (q(1)^2 - qv'*qv)*eye(3) + 2*(qv*qv') - 2*q(1)*skew_matrix(qv);
end
