function rec = estimator_scenario(p, seed, T, outage, noisy, bias, blind, opt)
% Truth and sensors for exercising src/estimation without Simulink: a rigid body with reaction
% wheels driven by a scripted wheel torque, disturbance torques, a gyro and a star tracker.
% The filter sees the commanded wheel torque and the wheel momentum; imperfect knowledge is set in
% opt (missing fields take the defaults below):
%   dist           scale of the disturbance torques (1)
%   torqueErr      relative error of the delivered wheel torque against the commanded one (0)
%   torqueBias     constant extra torque on the body that the filter is not told about, N m (0)
%   torqueNoise    white torque noise per step on the body, N m rms (0)
%   tachNoise      noise on the wheel momentum the filter receives, N m s rms (0)
%   inertiaErr     relative error of the filter's inertia diagonal against the truth (0)
%   rateDamp       rate feedback gain on the commanded wheel torque, 1/s (0.05): stands in for the
%                  attitude controller, which keeps the body rate bounded when a constant torque acts
%                  (set 0 for the pure scripted torque, where a constant torque spins the body up)
%   starRatioTrue    boresight / cross-boresight noise ratio of the simulated tracker (p value)
%   starRatioFilter  the same ratio as the filter assumes (p value)
% rec holds the estimation errors per step.
def = struct('dist', 1, 'torqueErr', 0, 'torqueBias', zeros(3, 1), 'torqueNoise', 0, ...
             'tachNoise', 0, 'inertiaErr', zeros(3, 1), 'starRatioTrue', [], 'starRatioFilter', [], 'rateDamp', 0.05);
names = fieldnames(def);
for k = 1:numel(names)
    if ~isfield(opt, names{k})
        opt.(names{k}) = def.(names{k});
    end
end
if isscalar(opt.inertiaErr)
    opt.inertiaErr = opt.inertiaErr*ones(3, 1);
end
pt = p;
if ~isempty(opt.starRatioTrue)
    pt.starBoresightRatio = opt.starRatioTrue;
end
starRoot = sqrtm(star_covariance(pt));
pf = p;
if ~isempty(opt.starRatioFilter)
    pf.starBoresightRatio = opt.starRatioFilter;
end
Dm = diag(1 + opt.inertiaErr(:));
pf.inertia = 0.5*(Dm*p.inertia + p.inertia*Dm);

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
wt0 = wt;
inertia_f = p.inertia;
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
        qm = quat_mul(qt, rotvec_to_quat(starRoot*randn(3, 1)));
        hMeas = hPrev;
        if opt.tachNoise > 0
            hMeas = hMeas + opt.tachNoise*randn(3, 1);
        end
    else
        wm = wt + b;
        qm = qt;
        hMeas = hPrev;
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
    x = mekf_step(x, wm, qm, new, tauPrev, hMeas, pf);
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
    tau = torque_profile(t) + opt.rateDamp*(inertia_f*(wt - wt0));
    hPrev = ht;
    extra = opt.torqueBias(:);
    if noisy && opt.torqueNoise > 0
        extra = extra + opt.torqueNoise*randn(3, 1);
    end
    [qt, wt, ht] = truth_step(qt, wt, ht, tau*(1 + opt.torqueErr), extra, t, dt, p, opt.dist);
    tauPrev = tau;
    if noisy
        b = b + randn(3, 1)*p.gyroRrw*sqrt(dt);
    end
end
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

function [q, w, h] = truth_step(q, w, h, tau, extra, t, dt, p, dist)
n = 2;
ds = dt/n;
for j = 1:n
    tj = t + (j - 1)*ds;
    d = dist*disturbance(tj + 0.5*ds) + extra;
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
