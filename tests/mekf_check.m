function mekf_check
% Checks of the multiplicative EKF in src/estimation without Simulink: exact propagation,
% covariance propagation against the exact matrix exponential, filter consistency and accuracy
% over Monte Carlo runs, and a 600 s star tracker outage. Sensors follow smallsat_params.
% Thresholds were fixed before the first run.
p = smallsat_params;
dt = p.gyroStep;
fails = 0;
total = 0;

const = @(t) [0.001; -0.0011; 0.0007];
rec = simulate(p, 0, 600, [], false, zeros(3, 1), const);
v = max(vecnorm(rec.e, 2, 1));
[fails, total] = report(fails, total, 'noise free attitude error', v, 1e-10, v < 1e-10);

w = [0.0012; -0.0011; 0.0016];
A = [-skew_matrix(w), -eye(3); zeros(3, 6)];
G = blkdiag(-eye(3), eye(3));
Qc = blkdiag(p.gyroArw^2*eye(3), p.gyroRrw^2*eye(3));
E = expm([-A, G*Qc*G'; zeros(6), A']*dt);
PhiEx = E(7:12, 7:12)';
QdEx = PhiEx*E(1:6, 7:12);
Aa = A*dt;
Phi = eye(6) + Aa + 0.5*(Aa*Aa);
su = p.gyroArw^2;
sv = p.gyroRrw^2;
Qd = [(su*dt + sv*dt^3/3)*eye(3), -(sv*dt^2/2)*eye(3); -(sv*dt^2/2)*eye(3), sv*dt*eye(3)];
v = max(abs(Phi - PhiEx), [], 'all');
[fails, total] = report(fails, total, 'Phi vs exact exponential', v, 1e-6, v < 1e-6);
v = max(abs(Qd - QdEx), [], 'all')/max(abs(QdEx), [], 'all');
[fails, total] = report(fails, total, 'Qd relative error vs Van Loan', v, 1e-3, v < 1e-3);

nRuns = 8;
nees = zeros(nRuns, 1);
rmsErr = zeros(nRuns, 1);
bz = zeros(nRuns, 3);
bErr = zeros(nRuns, 3);
for s = 1:nRuns
    rec = simulate(p, 100 + s - 1, 6000, [], true, p.gyroBias0, @truth_rate);
    sel = rec.t >= 600 & rec.newStar;
    idx = find(sel);
    ne = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        ne(k) = rec.e(:, idx(k))'*(rec.P3(:, :, idx(k)) \ rec.e(:, idx(k)));
    end
    nees(s) = mean(ne);
    rmsErr(s) = sqrt(mean(sum(rec.e(:, idx).^2, 1)/3));
    bErr(s, :) = rec.bErr(:, end)';
    bz(s, :) = rec.bErr(:, end)'./sqrt(diag(rec.Pb(:, :, end)))';
end
v = mean(nees);
[fails, total] = report(fails, total, 'NEES mean (ideal 3)', v, '2.0 .. 4.5', v > 2.0 && v < 4.5);
v = mean(rmsErr)/p.starSigma;
[fails, total] = report(fails, total, 'attitude rms / starSigma', v, '< 1', v < 1);
v = sqrt(mean(bz(:).^2));
[fails, total] = report(fails, total, 'bias error / filter sigma (rms)', v, '< 2', v < 2);
v = sqrt(mean(bErr(:).^2))/norm(p.gyroBias0);
[fails, total] = report(fails, total, 'bias error rms vs initial |bias|', v, '< 0.5', v < 0.5);

nOut = 4;
inside = zeros(nOut, 1);
worst = zeros(nOut, 1);
for s = 1:nOut
    rec = simulate(p, 200 + s - 1, 6000, [3000 3600], true, p.gyroBias0, @truth_rate);
    idx = find(rec.t >= 3000 & rec.t < 3600);
    ok = false(numel(idx), 1);
    for k = 1:numel(idx)
        ok(k) = all(abs(rec.e(:, idx(k))) < 3*sqrt(diag(rec.P3(:, :, idx(k)))));
    end
    inside(s) = mean(ok);
    worst(s) = max(vecnorm(rec.e(:, idx), 2, 1))/p.starSigma;
end
v = mean(inside);
[fails, total] = report(fails, total, 'outage: fraction inside 3 sigma', v, '> 0.97', v > 0.97);
fprintf('   worst outage attitude error / starSigma: %.2f\n', max(worst));
fprintf('   NEES per run: %s\n   rms/starSigma per run: %s\n', mat2str(nees', 3), mat2str(rmsErr'/p.starSigma, 3));
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
fprintf('%s  %-36s %10.3e (limit %s)\n', s, name, val, lim);
end

function w = truth_rate(t)
base = [0; -0.0622; 0]*pi/180;
osc = 0.02*pi/180*[sin(2*pi*t/700); cos(2*pi*t/450); sin(2*pi*t/300)];
slew = [0.3; 0; -0.2]*pi/180*(t >= 4500 && t < 4560);
w = base + osc + slew;
end

function rec = simulate(p, seed, T, outage, noisy, bias, rate)
rng(seed);
dt = p.gyroStep;
n = round(T/dt) + 1;
qt = [0.5; 0.5; -0.5; 0.5];
x = zeros(43, 1);
b = bias;
rec.t = zeros(1, n);
rec.e = zeros(3, n);
rec.P3 = zeros(3, 3, n);
rec.bErr = zeros(3, n);
rec.Pb = zeros(3, 3, n);
rec.newStar = false(1, n);
for k = 1:n
    t = (k - 1)*dt;
    w = rate(t);
    if noisy
        wm = w + b + randn(3, 1)*p.gyroArw/sqrt(dt);
        qm = quat_mul(qt, qexp(randn(3, 1)*p.starSigma));
    else
        wm = w + b;
        qm = qt;
    end
    new = mod(k - 1, round(p.starStep/dt)) == 0;
    if ~isempty(outage) && t >= outage(1) && t < outage(2)
        new = false;
    end
    if k == 1
        new = true;
    end
    x = mekf_step(x, wm, qm, new, p);
    d = quat_mul([x(1); -x(2:4)], qt);
    P = reshape(x(8:43), 6, 6);
    rec.t(k) = t;
    rec.e(:, k) = 2*d(2:4);
    rec.P3(:, :, k) = P(1:3, 1:3);
    rec.bErr(:, k) = b - x(5:7);
    rec.Pb(:, :, k) = P(4:6, 4:6);
    rec.newStar(k) = new;
    for j = 1:10
        qt = quat_mul(qt, qexp(rate(t + (j - 0.5)*dt/10)*dt/10));
        qt = qt/norm(qt);
    end
    if noisy
        b = b + randn(3, 1)*p.gyroRrw*sqrt(dt);
    end
end
end

function q = qexp(th)
a = norm(th);
if a > 1e-12
    q = [cos(a/2); sin(a/2)*th/a];
else
    q = [1; th/2];
end
end
