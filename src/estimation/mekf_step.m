function x = mekf_step(x, wMeas, qMeas, newStar, p)
% One gyro step of the multiplicative EKF, plus the star tracker update when newStar is nonzero.
% State x = [qHat(4); bHat(3); P(36)], error state [dtheta; db] with dtheta in the body frame
% (q = qHat (x) [1; dtheta/2]) and db = b - bHat. The filter starts itself from the first star
% tracker sample when x(1:4) is zero.
if norm(x(1:4)) == 0
    x = mekf_init(qMeas, p);
    return
end
dt = p.gyroStep;
q = x(1:4);
b = x(5:7);
P = reshape(x(8:43), 6, 6);

w = wMeas(:) - b;
nw = norm(w);
if nw*dt > 1e-12
    dq = [cos(0.5*nw*dt); sin(0.5*nw*dt)*w/nw];
else
    dq = [1; 0.5*dt*w];
end
q = quat_mul(q, dq);
q = q/norm(q);

A = [-skew_matrix(w), -eye(3); zeros(3, 6)]*dt;
Phi = eye(6) + A + 0.5*(A*A);
su = p.gyroArw^2;
sv = p.gyroRrw^2;
Qd = [(su*dt + sv*dt^3/3)*eye(3), -(sv*dt^2/2)*eye(3); ...
      -(sv*dt^2/2)*eye(3),         sv*dt*eye(3)];
P = Phi*P*Phi' + Qd;

if newStar ~= 0
    d = quat_mul([q(1); -q(2:4)], qMeas(:));
    if d(1) < 0
        d = -d;
    end
    z = 2*d(2:4);
    R = p.starSigma^2*eye(3);
    K = P(:, 1:3)/(P(1:3, 1:3) + R);
    dx = K*z;
    IKH = eye(6) - [K, zeros(6, 3)];
    P = IKH*P*IKH' + K*R*K';
    th = dx(1:3);
    nth = norm(th);
    if nth > 1e-12
        dq = [cos(0.5*nth); sin(0.5*nth)*th/nth];
    else
        dq = [1; 0.5*th];
    end
    q = quat_mul(q, dq);
    q = q/norm(q);
    b = b + dx(4:6);
end
P = 0.5*(P + P');
x = [q; b; P(:)];
end
