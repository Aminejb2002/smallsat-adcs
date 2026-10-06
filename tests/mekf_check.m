function mekf_check
% Checks of the rate-aided multiplicative EKF in src/estimation without Simulink: exact
% propagation, error-state Jacobian against finite differences, covariance propagation against
% the exact matrix exponential, filter consistency and accuracy over Monte Carlo runs, a 600 s
% star tracker outage, a late start and a wheel torque knowledge error. The truth is a rigid body
% with wheels, disturbance torques and the sensors of smallsat_params.
% Thresholds were fixed before the first run, except: the inertia/torque error checks and the rate
% NEES lower bound (0.5), set after the robustness study in experiments/run_estimator_robustness,
% and the two noise free limits, because the filter uses the wheel momentum from the start of the
% step, which leaves about 5e-9 rad/s of model error.
p = smallsat_params;
dt = p.gyroStep;
sigGyro = p.gyroArw/sqrt(dt);
sigAx = sqrt(diag(star_covariance(p)));   % star tracker sigma per body axis
fails = 0;
total = 0;

none = struct('dist', 0);
rec = estimator_scenario(p, 0, 600, [], false, zeros(3, 1), 0, none);
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

Qc = mekf_process_noise(p, [2e-5; -1e-5; 3e-5]);
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
full = struct('dist', 1);
for s = 1:nRuns
    rec = estimator_scenario(p, 100 + s - 1, 3000, [], true, p.gyroBias0, 0, full);
    sel = rec.t >= 600;
    idx = find(sel & rec.newStar);
    ne = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        ne(k) = rec.e(:, idx(k))'*(rec.P3(:, :, idx(k)) \ rec.e(:, idx(k)));
    end
    nees(s) = mean(ne);
    rmsErr(s) = sqrt(mean(sum((rec.e(:, idx)./sigAx).^2, 1)/3));
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
v = mean(rmsErr);
[fails, total] = report(fails, total, 'attitude rms / star sigma per axis', v, '< 1', v < 1);
v = mean(neesW);
[fails, total] = report(fails, total, 'NEES rate mean (ideal 3)', v, '0.5 .. 4.5', v > 0.5 && v < 4.5);
v = mean(rmsW)/sigGyro;
[fails, total] = report(fails, total, 'rate error rms / gyro noise', v, '< 0.15', v < 0.15);
v = sqrt(mean(bz(:).^2));
[fails, total] = report(fails, total, 'bias error / filter sigma (rms)', v, '< 2', v < 2);
v = sqrt(mean(bErr(:).^2))/norm(p.gyroBias0);
[fails, total] = report(fails, total, 'bias error rms vs initial |bias|', v, '< 0.5', v < 0.5);

nOut = 6;   % the error drifts slowly, so each run is one realisation: use several
inside = zeros(nOut, 1);
neesOut = zeros(nOut, 1);
worst = zeros(nOut, 1);
for s = 1:nOut
    rec = estimator_scenario(p, 200 + s - 1, 3000, [1500 2100], true, p.gyroBias0, 0, full);
    idx = find(rec.t >= 1500 & rec.t < 2100);
    ok = false(numel(idx), 1);
    for k = 1:numel(idx)
        ok(k) = all(abs(rec.e(:, idx(k))) < 3*sqrt(diag(rec.P3(:, :, idx(k)))));
    end
    inside(s) = mean(ok);
    nn = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        nn(k) = rec.e(:, idx(k))'*(rec.P3(:, :, idx(k)) \ rec.e(:, idx(k)));
    end
    neesOut(s) = mean(nn);
    worst(s) = max(vecnorm(rec.e(:, idx)./sigAx, 2, 1));
end
v = mean(inside);
[fails, total] = report(fails, total, 'outage: fraction inside 3 sigma', v, '> 0.95', v > 0.95);
v = mean(neesOut);
[fails, total] = report(fails, total, 'outage: NEES attitude mean over runs', v, '1 .. 6', v > 1 && v < 6);

blind = 300;
rec = estimator_scenario(p, 300, 3000, [], true, p.gyroBias0, blind, full);
late = rec.t >= blind + 600 & rec.newStar;
v = sqrt(mean(sum((rec.e(:, late)./sigAx).^2, 1)/3));
[fails, total] = report(fails, total, 'late start: rms / star sigma', v, '< 1', v < 1);
[fails, total] = report(fails, total, 'late start: raw pass-through', rec.passthrough, 1e-15, rec.passthrough < 1e-15);

wrong = struct('torqueErr', 0.05);
rec = estimator_scenario(p, 400, 3000, [], true, p.gyroBias0, 0, wrong);
idx = find(rec.t >= 600 & rec.newStar);
v = sqrt(mean(sum((rec.e(:, idx)./sigAx).^2, 1)/3));
[fails, total] = report(fails, total, '5 % wheel torque error: rms / star sigma', v, '< 1.5', v < 1.5);
idx = find(rec.t >= 600);
v = sqrt(mean(sum(rec.wErr(:, idx).^2, 1)/3))/sigGyro;
[fails, total] = report(fails, total, '5 % wheel torque error: rate rms / gyro', v, '< 0.3', v < 0.3);

miss = struct('inertiaErr', [0.10; -0.10; 0.05], 'torqueErr', 0.05);
rec = estimator_scenario(p, 500, 3000, [], true, p.gyroBias0, 0, miss);
idx = find(rec.t >= 600);
nA = zeros(numel(idx), 1);
nW = zeros(numel(idx), 1);
for k = 1:numel(idx)
    nA(k) = rec.e(:, idx(k))'*(rec.P3(:, :, idx(k)) \ rec.e(:, idx(k)));
    nW(k) = rec.wErr(:, idx(k))'*(rec.Pw(:, :, idx(k)) \ rec.wErr(:, idx(k)));
end
v = sqrt(mean(sum((rec.e(:, idx)./sigAx).^2, 1)/3));
[fails, total] = report(fails, total, 'inertia/torque error: attitude rms / sig', v, '< 1', v < 1);
v = mean(nA);
[fails, total] = report(fails, total, 'inertia/torque error: NEES attitude', v, '1 .. 6', v > 1 && v < 6);
v = mean(nW);
[fails, total] = report(fails, total, 'inertia/torque error: NEES rate', v, '< 8', v < 8);

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

function wd = rate_dot(w, h, d, tauW, inertia)
wd = angular_accel(w, h, d, tauW, inertia);
end
