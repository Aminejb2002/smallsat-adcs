function res = run_monte_carlo(nRuns, stopTime, only, inertiaSpread)
% Monte Carlo of the full mission (tumble, B-dot detumble, handover, wheel capture, pointing with
% magnetic momentum dumping) with the estimate in the loop. Per run: random tip-off rate (3 to
% 6 deg/s, random direction), random initial attitude, random initial gyro bias (+-1.5 deg/h per
% axis) and new noise seeds for magnetometer, gyro, gyro bias drift and star tracker.
% inertiaSpread > 0 also gives the plant its own inertia, I_true = D*I*D with D = diag(sqrt(1+e)), e uniform in
% +-spread per axis (so I_true stays symmetric and positive definite), drawn from a separate random stream so
% every other draw is the same as without the spread. The controller and filter keep the nominal inertia.
% Not varied: actuator and sensor error models, disturbance size. Draws use rng(2026), so the table is reproducible.
% Optional: stopTime (s) overrides the model stop time, only = run numbers to simulate.
% Statistics over the last 2000 s of each run. Logged every 10 s, so peaks between samples are not seen.
if nargin < 1, nRuns = 20; end
if nargin < 2, stopTime = []; end   % [] keeps the model's stop time
if nargin < 3 || isempty(only), only = 1:nRuns; end   % all draws are made anyway, so run i is the same case as in a full run
if nargin < 4, inertiaSpread = 0; end   % e.g. 0.10 for plus or minus 10 % per axis
root = fileparts(fileparts(mfilename('fullpath')));
p0 = smallsat_params;
period = 2*pi*sqrt((p0.Re + p0.altitude)^3/p0.mu);
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
cleanup = onCleanup(@() evalin('base', 'clear adcsCase paramOverride'));
assignin('base', 'adcsCase', 'tumble');
rng(2026);
rsI = RandStream('mt19937ar', 'Seed', 77);
seedBlocks = {sprintf('Random\nNumber'), [41 42 43]; sprintf('Random\nNumber1'), [51 52 53]; ...
              sprintf('Random\nNumber2'), [61 62 63]; 'noise', [11 22 33]};
arc = 180/pi*3600;
res = struct('run', num2cell(1:nRuns));
fprintf('%4s %8s %9s %9s %9s %9s %9s %9s %9s %7s %7s\n', 'run', 'tip-off', 'handover', 'capture', 'point', 'point', 'knowledge', 'h last', 'dipole', 'NEES', 'dI');
fprintf('%4s %8s %9s %9s %9s %9s %9s %9s %9s %7s %7s\n', '', '[deg/s]', '[s]', '[s]', 'mean[deg]', 'max[deg]', 'rms[arcs]', 'max[Nms]', 'max[Am2]', 'mean', 'max[%]');
for i = 1:nRuns
    d = randn(3, 1); d = d/norm(d);
    rate = (3 + 3*rand)*pi/180;
    qi = randn(4, 1); qi = qi/norm(qi);
    bias = (2*rand(3, 1) - 1)*1.5*pi/180/3600;
    eI = (2*rand(rsI, 3, 1) - 1)*inertiaSpread;   % drawn for every run, also the skipped ones
    if ~ismember(i, only), continue; end
    ov = struct('tipOffRate', rate*d, 'initialQuat', qi, 'gyroBias0', bias);
    if inertiaSpread > 0
        D = diag(sqrt(1 + eI));
        ov.inertiaTruth = D*p0.inertia*D;
    end
    assignin('base', 'paramOverride', ov);
    in = Simulink.SimulationInput(mdl);
    in = in.setBlockParameter([mdl '/use_estimate'], 'Value', '1');
    if ~isempty(stopTime), in = in.setModelParameter('StopTime', num2str(stopTime)); end
    for b = 1:size(seedBlocks, 1)
        in = in.setBlockParameter([mdl '/Sensors/' seedBlocks{b, 1}], 'Seed', mat2str(seedBlocks{b, 2} + 1000*i));
    end
    out = sim(in);
    m = metrics(out, p0, period, arc);
    m.tipOff = rate*180/pi;
    m.inertiaErr = max(abs(eI))*100;
    fn = fieldnames(m);
    for k = 1:numel(fn), res(i).(fn{k}) = m.(fn{k}); end
    fprintf('%4d %8.2f %9.0f %9.0f %9.4f %9.4f %9.1f %9.3f %9.1f %7.2f %7.1f\n', i, m.tipOff, m.handover, m.capture, ...
        m.pointMean, m.pointMax, m.knowRms, m.hMax, m.dipoleMax, m.nees, m.inertiaErr);
end
res = res(only);
ok = [res.handover] > 0 & ~isnan([res.capture]);
fprintf('\n%d of %d runs reached pointing (handover and capture inside the run)\n', sum(ok), nRuns);
f = {'handover', 'capture', 'pointMean', 'pointMax', 'knowRms', 'hMax', 'dipoleMax', 'nees', 'inside3'};
fprintf('%-10s %12s %12s %12s %12s\n', '', 'mean', 'std', 'worst', 'best');
for k = 1:numel(f)
    v = [res(ok).(f{k})];
    fprintf('%-10s %12.4g %12.4g %12.4g %12.4g\n', f{k}, mean(v), std(v), max(v), min(v));
end
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
if inertiaSpread > 0, name = 'monte_carlo_inertia.mat';
elseif isempty(stopTime) && numel(only) == nRuns, name = 'monte_carlo.mat'; else, name = 'monte_carlo_rerun.mat'; end
save(fullfile(outDir, name), 'res');
if ~wasLoaded, close_system(mdl, 0); end
end

function m = metrics(out, p, period, arc)
t = out.t_log(:);
q = log_rows(out.q_log);
qref = log_rows(out.qref_log);
tauW = log_rows(out.tauW_log);
h = log_rows(out.h_log);
mt = log_rows(out.m_log);
errDeg = zeros(numel(t), 1);
for k = 1:numel(t)
    qe = quat_mul([qref(k, 1); -qref(k, 2:4)'], q(k, :)');
    errDeg(k) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
end
i0 = find(any(tauW ~= 0, 2), 1);
if isempty(i0)
    m = struct('handover', 0, 'capture', NaN, 'pointMean', NaN, 'pointMax', NaN, 'knowRms', NaN, 'hMax', NaN, 'dipoleMax', NaN, 'nees', NaN, 'inside3', NaN);
    return
end
tEng = t(i0);
bad = find(errDeg >= 0.1 & t >= tEng, 1, 'last');
if isempty(bad), cap = 0; elseif bad == numel(t), cap = NaN; else, cap = t(bad + 1) - tEng; end
tail = t >= t(end) - 2000;
qHat = log_rows(out.estQ_log);
step = round(10/p.gyroStep);
idx = 1:step:size(qHat, 1);
qHat = qHat(idx, :);
sig = log_rows(out.estSig_log);
sig = sig(idx, :);
e = zeros(sum(tail), 1);
ev = zeros(sum(tail), 3);
rows = find(tail);
for j = 1:numel(rows)
    d = quat_mul([qHat(rows(j), 1); -qHat(rows(j), 2:4)'], q(rows(j), :)');
    e(j) = 2*atan2(norm(d(2:4)), abs(d(1)))*arc;
    if d(1) < 0, d = -d; end
    ev(j, :) = 2*d(2:4)';
end
lastOrbit = t >= t(end) - period;
after = t >= tEng + 10;
m.handover = tEng;
m.capture = cap;
m.pointMean = mean(errDeg(tail));
m.pointMax = max(errDeg(tail));
m.knowRms = sqrt(mean(e.^2));
m.nees = mean(sum((ev./sig(rows, :)).^2, 2));       % ideal 3
m.inside3 = mean(all(abs(ev) < 3*sig(rows, :), 2));
m.hMax = max(vecnorm(h(lastOrbit, :), 2, 2));
m.dipoleMax = max(max(abs(mt(after, :))));
end
