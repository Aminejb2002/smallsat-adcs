function reference_check
% Checks of the reference models: quaternion math, J2 orbit, rigid body,
% gravity-gradient libration and the disturbance budget against independent values.

p = smallsat_params;
tests = {
    'quaternion round trip and orthonormality', @() test_roundtrip()
    'quaternion rate matches DCM rate',         @() test_kinematics()
    'quat_to_dcm matches Aerospace quat2dcm',   @() test_toolbox_convention()
    'J2 acceleration is the potential gradient', @() test_potential_gradient(p)
    'Sun-synchronous inclination at 600 km',    @() test_sso(p)
    'J2 energy conservation, 12 h',             @() test_energy(p)
    'nodal precession equals SSO rate, 5 days', @() test_nodal(p)
    'torque-free momentum and energy',          @() test_torque_free()
    'axisymmetric body rate precession',        @() test_axisymmetric()
    'gravity-gradient pitch libration',         @() test_libration(p)
    'disturbance budget vs independent values', @() test_budget()
    'wheel torque conserves total momentum',    @() test_wheel_momentum()
    'nadir reference matches nadir_dcm',        @() test_nadir_reference(p)
    'pointing law is zero at zero error',       @() test_pointing_zero(p)
    'pointing law reduces to PD for small error', @() test_pointing_linear(p)
    };

fails = 0;
for k = 1:size(tests, 1)
    [ok, detail] = tests{k, 2}();
    if ok
        status = 'PASS';
    else
        status = 'FAIL';
        fails = fails + 1;
    end
    fprintf('%s  %-44s %s\n', status, tests{k, 1}, detail);
end
fprintf('\n%d of %d checks passed\n', size(tests, 1) - fails, size(tests, 1));
if fails > 0
    error('reference_check:failed', '%d checks failed', fails);
end
end

function [ok, detail] = test_roundtrip
rng(1);
worst = 0;
for k = 1:200
    q = randn(4, 1);
    q = q/norm(q);
    A = quat_to_dcm(q);
    if q(1) < 0, q = -q; end
    worst = max([worst, max(abs(A*A' - eye(3)), [], 'all'), abs(det(A) - 1), ...
        max(abs(dcm_to_quat(A) - q))]);
end
ok = worst < 1e-10;
detail = sprintf('worst error %.2e', worst);
end

function [ok, detail] = test_kinematics
rng(2);
q = randn(4, 1);
q = q/norm(q);
w = [0.03; -0.05; 0.02];
dt = 1e-6;
dA = (quat_to_dcm(q + dt*quat_rate(q, w)) - quat_to_dcm(q - dt*quat_rate(q, w)))/(2*dt);
err = max(abs(dA + skew_matrix(w)*quat_to_dcm(q)), [], 'all');
ok = err < 1e-8;
detail = sprintf('error %.2e', err);
end

function [ok, detail] = test_toolbox_convention
rng(4);
worst = 0;
for k = 1:50
    q = randn(1, 4);
    q = q/norm(q);
    worst = max(worst, max(abs(quat_to_dcm(q) - quat2dcm(q)), [], 'all'));
end
ok = worst < 1e-12;
detail = sprintf('worst difference %.2e', worst);
end

function [ok, detail] = test_potential_gradient(p)
r = [5.2e6; -3.9e6; 3.1e6];
grad = zeros(3, 1);
for i = 1:3
    d = zeros(3, 1);
    d(i) = 1;
    grad(i) = (j2_potential(r + d, p) - j2_potential(r - d, p))/2;
end
acc = j2_accel(r, p);
err = max(abs(acc + grad))/max(abs(acc));
ok = err < 1e-6;
detail = sprintf('relative error %.2e', err);
end

function [ok, detail] = test_sso(p)
incDeg = sso_inclination(p.Re + 600e3, p)*180/pi;
ok = incDeg > 97.78 && incDeg < 97.80;
detail = sprintf('%.3f deg', incDeg);
end

function [t, x] = propagate(a, inc, raan, tspan, p, tol)
[r0, v0] = circular_state(a, inc, raan, 0, p);
opts = odeset('RelTol', tol, 'AbsTol', 1e-6);
[t, x] = ode113(@(t, x) [x(4:6); j2_accel(x(1:3), p)], tspan, [r0; v0], opts);
end

function [ok, detail] = test_energy(p)
a = p.Re + 600e3;
[~, x] = propagate(a, sso_inclination(a, p), 0.3, 0:30:43200, p, 1e-12);
e = zeros(size(x, 1), 1);
for k = 1:size(x, 1)
    e(k) = 0.5*(x(k, 4:6)*x(k, 4:6)') + j2_potential(x(k, 1:3)', p);
end
err = max(abs(e - e(1)))/abs(e(1));
ok = err < 1e-9;
detail = sprintf('relative drift %.2e', err);
end

function [ok, detail] = test_nodal(p)
a = p.Re + 600e3;
[t, x] = propagate(a, sso_inclination(a, p), 0.3, 0:60:5*86400, p, 1e-11);
h = cross(x(:, 1:3), x(:, 4:6), 2);
raan = unwrap(atan2(h(:, 1), -h(:, 2)));
c = polyfit(t, raan, 1);
target = 2*pi/p.yearSeconds;
err = abs(c(1) - target)/target;
ok = err < 0.01;
detail = sprintf('%.4f deg/day, error %.2f %%', c(1)*86400*180/pi, 100*err);
end

function [ok, detail] = test_torque_free
inertia = [26 0.6 -0.4; 0.6 24 0.5; -0.4 0.5 16];
rng(3);
q0 = randn(4, 1);
q0 = q0/norm(q0);
x0 = [q0; 0.05; -0.08; 0.12];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, x] = ode113(@(t, x) rigid_body_rhs(x, inertia, zeros(3, 1)), linspace(0, 200, 401), x0, opts);
n = size(x, 1);
H = zeros(n, 3);
ke = zeros(n, 1);
for k = 1:n
    w = x(k, 5:7)';
    H(k, :) = (quat_to_dcm(x(k, 1:4))'*(inertia*w))';
    ke(k) = 0.5*w'*inertia*w;
end
errH = max(vecnorm(H - H(1, :), 2, 2))/norm(H(1, :));
errE = max(abs(ke - ke(1)))/ke(1);
errQ = max(abs(vecnorm(x(:, 1:4), 2, 2) - 1));
ok = errH < 1e-9 && errE < 1e-9 && errQ < 1e-9;
detail = sprintf('H %.1e, energy %.1e, |q| %.1e', errH, errE, errQ);
end

function [ok, detail] = test_axisymmetric
inertia = diag([20 20 30]);
wz = 0.5;
wxy = 0.1;
x0 = [1; 0; 0; 0; wxy; 0; wz];
t = linspace(0, 40, 201);
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, x] = ode113(@(t, x) rigid_body_rhs(x, inertia, zeros(3, 1)), t, x0, opts);
lam = (30 - 20)/20*wz;
err = max([max(abs(x(:, 5) - wxy*cos(lam*t(:)))), max(abs(x(:, 6) - wxy*sin(lam*t(:)))), ...
    max(abs(x(:, 7) - wz))]);
ok = err < 1e-8;
detail = sprintf('error %.2e', err);
end

function r = orbit_r(t, a, inc, n, p)
r = circular_state(a, inc, 0.3, n*t, p);
end

function [ok, detail] = test_libration(p)
a = p.Re + 600e3;
n = sqrt(p.mu/a^3);
inertia = diag([26 24 16]);
inc = 97.8*pi/180;
[r0, v0] = circular_state(a, inc, 0.3, 0, p);
th0 = 2*pi/180;
ry = [cos(th0) 0 -sin(th0); 0 1 0; sin(th0) 0 cos(th0)];
A0 = ry*nadir_dcm(r0, v0);
hHat = cross(r0, v0)/norm(cross(r0, v0));
x0 = [dcm_to_quat(A0); A0*(n*hHat)];

torque = @(t, x) torque_gravity_gradient(quat_to_dcm(x(1:4))*orbit_r(t, a, inc, n, p), inertia, a, p);
tEnd = 4*2*pi/n;
opts = odeset('RelTol', 1e-11, 'AbsTol', 1e-13);
[t, x] = ode113(@(t, x) rigid_body_rhs(x, inertia, torque(t, x)), 0:5:tEnd, x0, opts);

pitch = zeros(numel(t), 1);
for k = 1:numel(t)
    [r, v] = circular_state(a, inc, 0.3, n*t(k), p);
    Aerr = quat_to_dcm(x(k, 1:4))*nadir_dcm(r, v)';
    pitch(k) = Aerr(1, 3);
end
s = sign(pitch);
idx = find(s(1:end-1) ~= s(2:end));
period = 2*mean(diff(t(idx)));
expected = 2*pi/(n*sqrt(3*(26 - 16)/24));
err = abs(period - expected)/expected;
ok = err < 0.005 && max(abs(pitch)) < sin(th0)*1.05;
detail = sprintf('period %.1f s vs %.1f s (%.2f %%)', period, expected, 100*err);
end

function [ok, detail] = test_budget
b = disturbance_budget;
expected = [2.253e-6, 8.606e-7, 4.398e-7, 1.783e-5];
err = max(abs(b.peak(1:4) - expected)./expected);
ok = err < 2e-3 && abs(b.inclinationDeg - 97.788) < 2e-3 && abs(b.periodSeconds - 5801.2) < 0.1;
detail = sprintf('peak torques within %.2f %% of reference', 100*err);
end

function [ok, detail] = test_wheel_momentum
inertia = [26 0.6 -0.4; 0.6 24 0.5; -0.4 0.5 16];
rng(5);
q0 = randn(4, 1);
q0 = q0/norm(q0);
x0 = [q0; 0.02; -0.01; 0.03; 0.1; -0.2; 0.05];
tauW = @(t) 0.01*[sin(0.1*t); cos(0.07*t); sin(0.03*t)];
f = @(t, x) [quat_rate(x(1:4), x(5:7)); ...
    angular_accel(x(5:7), x(8:10), zeros(3, 1), tauW(t), inertia); tauW(t)];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, x] = ode113(f, linspace(0, 300, 301), x0, opts);
H = zeros(size(x, 1), 3);
for k = 1:size(x, 1)
    H(k, :) = (quat_to_dcm(x(k, 1:4))'*(inertia*x(k, 5:7)' + x(k, 8:10)'))';
end
err = max(vecnorm(H - H(1, :), 2, 2))/norm(H(1, :));
ok = err < 1e-9;
detail = sprintf('relative drift %.2e', err);
end

function [ok, detail] = test_nadir_reference(p)
a = p.Re + p.altitude;
[r, v] = circular_state(a, sso_inclination(a, p), 0.3, 1.1, p);
[qref, wref] = nadir_reference(r, v);
errA = max(abs(quat_to_dcm(qref) - nadir_dcm(r, v)), [], 'all');
errW = abs(wref(2) + sqrt(p.mu/a^3))/sqrt(p.mu/a^3);
ok = errA < 1e-12 && errW < 1e-12 && wref(1) == 0 && wref(3) == 0;
detail = sprintf('DCM %.1e, rate %.1e', errA, errW);
end

function [ok, detail] = test_pointing_zero(p)
a = p.Re + p.altitude;
[r, v] = circular_state(a, sso_inclination(a, p), 0.3, 0.7, p);
[qref, wref] = nadir_reference(r, v);
u = pointing_law(qref, wref, qref, wref, p);
ok = norm(u) < 1e-12;
detail = sprintf('|u| %.1e N m', norm(u));
end

function [ok, detail] = test_pointing_linear(p)
a = p.Re + p.altitude;
[r, v] = circular_state(a, sso_inclination(a, p), 0.3, 0.7, p);
[qref, wref] = nadir_reference(r, v);
th = 1*pi/180;
rx = [1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];
q = dcm_to_quat(rx*nadir_dcm(r, v));
Ae = quat_to_dcm(q)*nadir_dcm(r, v)';
u = pointing_law(q, Ae*wref, qref, wref, p);
s = size_pointing(p);
err = abs(abs(u(1)) - s.kp(1)*th)/(s.kp(1)*th);
ok = err < 1e-3 && abs(u(2)) < 1e-6*abs(u(1)) && abs(u(3)) < 1e-6*abs(u(1));
detail = sprintf('torque %.3e vs kp*theta %.3e', abs(u(1)), s.kp(1)*th);
end
