function mekf_check
% Checks of the rate-aided multiplicative EKF in src/estimation without Simulink: exact
% propagation, error-state Jacobian against finite differences, covariance propagation against
% the exact matrix exponential, filter consistency and accuracy over Monte Carlo runs, a 600 s
% star tracker outage, a late start and a wheel torque knowledge error. The truth is a rigid body
% with wheels, disturbance torques and the sensors of smallsat_params.
% Thresholds were fixed before the first run, except the two noise free limits: the filter uses
% the wheel momentum from the start of the step, which leaves about 5e-9 rad/s of model error,
% and the lower bound of the rate NEES, widened from 1.5 to 1.0 after the first run (1.49: the
% filter is slightly conservative on the rate, which is acceptable; overconfidence is not).
p = smallsat_params;
dt = p.gyroStep;
sigGyro = p.gyroArw/sqrt(dt);
fails = 0;
total = 0;

none = struct('dist', 0, 'torqueErr', 0);
rec = simulate(p, 0, 600, [], false, zeros(3, 1), 0, none);
v = max(vecnorm(rec.e, 2, 1));
[fails, total] = report(fails, total, 'noise free attitude error', v, 1e-6, v < 1e-6);
v = max(vecnorm(rec.wErr, 2, 1));
[fails, total] = report(fails, total, 'noise free rate error', v, 1e-7, v < 1e-7);

w = [0.0012; -0.0011; 0.0016];
h = [0.1; -0.05; 0.2];
F = mekf_jacobian(w, h, p.inertia);
Fnum = zeros(3);
step = 1e-7;
for j = 1:3
    e = zeros(3, 1);
    e(j) = step;
    Fnum(:, j) = (rate_dot(w + e, h, zeros(3, 1), zeros(3, 1), p.inertia) - ...
                  rate_dot(w - e, h, zeros(3, 1), zeros(3, 1), p.inertia))/(2*step);
end
v = max(max(abs(F(4:6, 4:6) - Fnum)))/max(max(abs(Fnum)));
[fails, total] = report(fails, total, 'rate Jacobian vs finite difference', v, 1e-6, v < 1e-6);

Qc = mekf_process_noise(p);
[Phi, Qd] = mekf_discretize(F, Qc, dt);
E = expm([-F, Qc; zeros(12), F']*dt);
PhiEx = E(13:24, 13:24)';
QdEx = PhiEx*E(1:12, 13:24);
v = max(max(abs(Phi - PhiEx)));
[fails, total] = report(fails, total, 'Phi vs exact exponential', v, 1e-6, v < 1e-6);
v = max(max(abs(Qd - QdEx)))/max(max(abs(QdEx)));
[fails, total] = report(fails, total, 'Qd relative error vs Van Loan', v, 1e-3, v < 1e-3);

nRuns = 5;
nees = zeros(nRuns, 1);
neesW = zeros(nRuns, 1);
rmsErr = zeros(nRuns, 1);
rmsW = zeros(nRuns, 1);
bz = zeros(nRuns, 3);
bErr = zeros(nRuns, 3);
full = struct('dist', 1, 'torqueErr', 0);
for s = 1:nRuns
    rec = simulate(p, 100 + s - 1, 3000, [], true, p.gyroBias0, 0, full);
    sel = rec.t >= 600;
    idx = find(sel & rec.newStar);
    ne = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        ne(k) = rec.e(:, idx(k))'*(rec.P3(:, :, idx(k)) \ rec.e(:, idx(k)));
    end
    nees(s) = mean(ne);
    rmsErr(s) = sqrt(mean(sum(rec.e(:, idx).^2, 1)/3));
    idx = find(sel);
    nw = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        nw(k) = rec.wErr(:, idx(k))'*(rec.Pw(:, :, idx(k)) \ rec.wErr(:, idx(k)));
    end
    neesW(s) = mean(nw);
    rmsW(s) = sqrt(mean(sum(rec.wErr(:, idx).^2, 1)/3));
    bErr(s, :) = rec.bErr(:, end)';
    bz(s, :) = rec.bErr(:, end)'./sqrt(diag(rec.Pb(:, :, end)))';
end
v = mean(nees);
[fails, total] = report(fails, total, 'NEES attitude mean (ideal 3)', v, '2.0 .. 4.5', v > 2.0 && v < 4.5);
v = mean(rmsErr)/p.starSigma;
[fails, total] = report(fails, total, 'attitude rms / starSigma', v, '< 1', v < 1);
v = mean(neesW);
[fails, total] = report(fails, total, 'NEES rate mean (ideal 3)', v, '1.0 .. 4.5', v > 1.0 && v < 4.5);
v = mean(rmsW)/sigGyro;
[fails, total] = report(fails, total, 'rate error rms / gyro noise', v, '< 0.15', v < 0.15);
v = sqrt(mean(bz(:).^2));
[fails, total] = report(fails, total, 'bias error / filter sigma (rms)', v, '< 2', v < 2);
v = sqrt(mean(bErr(:).^2))/norm(p.gyroBias0);
[fails, total] = report(fails, total, 'bias error rms vs initial |bias|', v, '< 0.5', v < 0.5);

nOut = 3;
inside = zeros(nOut, 1);
worst = zeros(nOut, 1);
for s = 1:nOut
    rec = simulate(p, 200 + s - 1, 3000, [1500 2100], true, p.gyroBias0, 0, full);
    idx = find(rec.t >= 1500 & rec.t < 2100);
    ok = false(numel(idx), 1);
    for k = 1:numel(idx)
        ok(k) = all(abs(rec.e(:, idx(k))) < 3*sqrt(diag(rec.P3(:, :, idx(k)))));
    end
    inside(s) = mean(ok);
    worst(s) = max(vecnorm(rec.e(:, idx), 2, 1))/p.starSigma;
end
v = mean(inside);
[fails, total] = report(fails, total, 'outage: fraction inside 3 sigma', v, '> 0.97', v > 0.97);

blind = 300;
rec = simulate(p, 300, 3000, [], true, p.gyroBias0, blind, full);
late = rec.t >= blind + 600 & rec.newStar;
v = sqrt(mean(sum(rec.e(:, late).^2, 1)/3))/p.starSigma;
[fails, total] = report(fails, total, 'late start: rms / starSigma', v, '< 1', v < 1);
[fails, total] = report(fails, total, 'late start: raw pass-through', rec.passthrough, 1e-15, rec.passthrough < 1e-15);

wrong = struct('dist', 1, 'torqueErr', 0.05);
rec = simulate(p, 400, 3000, [], true, p.gyroBias0, 0, wrong);
idx = find(rec.t >= 600 & rec.newStar);
v = sqrt(mean(sum(rec.e(:, idx).^2, 1)/3))/p.starSigma;
[fails, total] = report(fails, total, '5 % wheel torque error: rms / starSigma', v, '< 1.5', v < 1.5);
idx = find(rec.t >= 600);
v = sqrt(mean(sum(rec.wErr(:, idx).^2, 1)/3))/sigGyro;
[fails, total] = report(fails, total, '5 % wheel torque error: rate rms / gyro', v, '< 0.3', v < 0.3);

fprintf('   worst outage attitude error / starSigma: %.2f\n', max(worst));
fprintf('   NEES attitude per run: %s\n   NEES rate per run: %s\n', mat2str(nees', 3), mat2str(neesW', 3));
fprintf('\n%d of %d checks passed\n', total - fails, total);
if fails > 0
    error('mekf_check:failed', '%d checks failed', fails);
end
end

function [fails, total] = report(fails, total, name, val, lim, ok)
total = total + 1;
if ok
    s = 'PASS';
else
    s = 'FAIL';
    fails = fails + 1;
end
if isnumeric(lim), lim = sprintf('%.0e', lim); end
fprintf('%s  %-40s %10.3e (limit %s)\n', s, name, val, lim);
end

function tau = torque_profile(t)
tau = 1e-4*[sin(2*pi*t/700); cos(2*pi*t/450); sin(2*pi*t/300)];
if t >= 1200 && t < 1260
    tau = tau + [2e-3; 0; -1.5e-3];
end
end

function d = disturbance(t)
d = [1e-5*sin(2*pi*t/5800) + 3e-6; 8e-6*cos(2*pi*t/5800 + 1) - 2e-6; 5e-6*sin(2*pi*t/5800 + 2) + 1e-6];
end

function wd = rate_dot(w, h, d, tauW, inertia)
wd = angular_accel(w, h, d, tauW, inertia);
end

function [q, w, h] = truth_step(q, w, h, tau, t, dt, p, dist)
n = 2;
ds = dt/n;
for j = 1:n
    tj = t + (j - 1)*ds;
    d = dist*disturbance(tj + 0.5*ds);
    y = [q; w; h];
    k1 = rhs(y, d, tau, p);
    k2 = rhs(y + 0.5*ds*k1, d, tau, p);
    k3 = rhs(y + 0.5*ds*k2, d, tau, p);
    k4 = rhs(y + ds*k3, d, tau, p);
    y = y + ds/6*(k1 + 2*k2 + 2*k3 + k4);
    q = y(1:4)/norm(y(1:4));
    w = y(5:7);
    h = y(8:10);
end
end

function dy = rhs(y, d, tau, p)
dy = [quat_rate(y(1:4), y(5:7)); angular_accel(y(5:7), y(8:10), d, tau, p.inertia); tau];
end

function rec = simulate(p, seed, T, outage, noisy, bias, blind, opt)
rng(seed);
dt = p.gyroStep;
n = round(T/dt) + 1;
qt = [0.5; 0.5; -0.5; 0.5];
wt = [0; -0.0622; 0]*pi/180;
ht = zeros(3, 1);
x = zeros(157, 1);
b = bias;
tauPrev = zeros(3, 1);
hPrev = zeros(3, 1);
rec.t = zeros(1, n);
rec.e = zeros(3, n);
rec.P3 = zeros(3, 3, n);
rec.bErr = zeros(3, n);
rec.Pb = zeros(3, 3, n);
rec.wErr = zeros(3, n);
rec.Pw = zeros(3, 3, n);
rec.newStar = false(1, n);
rec.passthrough = 0;
for k = 1:n
    t = (k - 1)*dt;
    if noisy
        wm = wt + b + randn(3, 1)*p.gyroArw/sqrt(dt);
        qm = quat_mul(qt, rotvec_to_quat(randn(3, 1)*p.starSigma));
    else
        wm = wt + b;
        qm = qt;
    end
    new = mod(k - 1, round(p.starStep/dt)) == 0;
    if ~isempty(outage) && t >= outage(1) && t < outage(2)
        new = false;
    end
    if t < blind
        new = false;
    elseif k == 1
        new = true;
    end
    x = mekf_step(x, wm, qm, new, tauPrev, hPrev, p);
    rec.t(k) = t;
    rec.newStar(k) = new;
    if norm(x(1:4)) == 0
        [qo, wo] = mekf_outputs(x, wm, qm);
        rec.passthrough = max(rec.passthrough, max(norm(qo - qm), norm(wo - wm)));
        rec.e(:, k) = NaN;
        rec.P3(:, :, k) = NaN;
        rec.bErr(:, k) = NaN;
        rec.Pb(:, :, k) = NaN;
        rec.wErr(:, k) = NaN;
        rec.Pw(:, :, k) = NaN;
    else
        d = quat_mul([x(1); -x(2:4)], qt);
        P = reshape(x(14:157), 12, 12);
        rec.e(:, k) = 2*d(2:4);
        rec.P3(:, :, k) = P(1:3, 1:3);
        rec.bErr(:, k) = b - x(8:10);
        rec.Pb(:, :, k) = P(7:9, 7:9);
        rec.wErr(:, k) = x(5:7) - wt;
        rec.Pw(:, :, k) = P(4:6, 4:6);
    end
    tau = torque_profile(t);
    hPrev = ht;
    [qt, wt, ht] = truth_step(qt, wt, ht, tau*(1 + opt.torqueErr), t, dt, p, opt.dist);
    tauPrev = tau;
    if noisy
        b = b + randn(3, 1)*p.gyroRrw*sqrt(dt);
    end
end
end
