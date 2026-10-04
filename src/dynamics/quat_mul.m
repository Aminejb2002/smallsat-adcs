function r = quat_mul(p, q)
p = p(:);
q = q(:);
r = [p(1)*q(1) - p(2:4)'*q(2:4); ...
     p(1)*q(2:4) + q(1)*p(2:4) + cross(p(2:4), q(2:4))];
end
