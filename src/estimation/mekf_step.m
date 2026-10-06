function x = mekf_step(x, wMeas, qMeas, newStar, tauW, h, p)
% One gyro step of the multiplicative EKF with a rate and disturbance torque model.
% State x = [qHat(4); wHat(3); bHat(3); dHat(3); P(144)], error state [dtheta; dw; db; dd] with
% dtheta in the body frame (q = qHat (x) [1; dtheta/2]). The rate is propagated with
%   inertia*wdot = d - tauW - w x (inertia*w + h)
% from the wheel torque tauW applied over the step and the wheel momentum h, and corrected by
% the gyro (wMeas = w + b + noise) every step and by the star tracker when newStar is nonzero.
% The filter starts from the first star fix after the state was zero and stays at zero until then.
if norm(x(1:4)) == 0
    if newStar ~= 0
        x = mekf_init(qMeas, wMeas, p);
    end
    return
end
dt = p.gyroStep;
inertia = p.inertia;
tauW = tauW(:);
h = h(:);
q = x(1:4);
w = x(5:7);
b = x(8:10);
d = x(11:13);
P = reshape(x(14:157), 12, 12);

wNew = w + dt*(inertia\(d - tauW - cross(w, inertia*w + h)));
wBar = 0.5*(w + wNew);
q = quat_mul(q, rotvec_to_quat(wBar*dt));
q = q/norm(q);
[Phi, Qd] = mekf_discretize(mekf_jacobian(wBar, h, inertia), mekf_process_noise(p), dt);
P = Phi*P*Phi' + Qd;
w = wNew;

Hg = [zeros(3), eye(3), eye(3), zeros(3)];
[q, w, b, d, P] = kalman_update(q, w, b, d, P, Hg, (p.gyroArw^2/dt)*eye(3), wMeas(:) - w - b);
if newStar ~= 0
    e = quat_mul([q(1); -q(2:4)], qMeas(:));
    if e(1) < 0
        e = -e;
    end
    Hs = [eye(3), zeros(3, 9)];
    [q, w, b, d, P] = kalman_update(q, w, b, d, P, Hs, p.starSigma^2*eye(3), 2*e(2:4));
end
P = 0.5*(P + P');
x = [q; w; b; d; P(:)];
end

function [q, w, b, d, P] = kalman_update(q, w, b, d, P, H, R, z)
S = H*P*H' + R;
K = (P*H')/S;
dx = K*z;
IKH = eye(12) - K*H;
P = IKH*P*IKH' + K*R*K';
q = quat_mul(q, rotvec_to_quat(dx(1:3)));
q = q/norm(q);
w = w + dx(4:6);
b = b + dx(7:9);
d = d + dx(10:12);
end
