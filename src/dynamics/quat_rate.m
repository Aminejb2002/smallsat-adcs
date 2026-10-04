function dq = quat_rate(q, w)
dq = 0.5*quat_mul(q, [0; w(:)]);
end
