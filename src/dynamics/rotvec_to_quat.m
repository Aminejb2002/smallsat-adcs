function q = rotvec_to_quat(th)
% Scalar-first quaternion of the rotation vector th (rad).
th = th(:);
n = norm(th);
if n > 1e-12
    q = [cos(0.5*n); sin(0.5*n)*th/n];
else
    q = [1; 0.5*th];
end
end
