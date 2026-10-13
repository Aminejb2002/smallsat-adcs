function res = run_detumble_interference_sim(nRuns, c, nCycle, nOff, only, duration, kScale, tauS)
% The detumble interference study in the Simulink copy (smallsat_adcs_interference.slx). Same arguments and
% same initial conditions as experiments/run_detumble_interference.m (plain MATLAB), so the two can be compared:
% run i has the tip-off rate, direction and attitude of run i in run_monte_carlo.
%   c        coupling torquer dipole -> magnetometer, T per A m^2
%   nCycle   B-dot cycle, s;  nOff  torquers off at the start of each cycle, s
%   kScale   factor on the B-dot gain;  tauS  settling time constant of the coil field, s
% Detumble time = hand-over to the wheels (first non-zero wheel torque), as in run_monte_carlo; 10 s log.
% After the hand-over the slew to nadir raises the body rate again, so 'last time above 0.5 deg/s' would be wrong here.
% Peak dipole: before the hand-over only.
% Only the B-dot phase has the cycle and the leakage; momentum dumping and pointing are as in the original.
if nargin < 1 || isempty(nRuns), nRuns = 30; end
if nargin < 2 || isempty(c), c = 0; end
if nargin < 3 || isempty(nCycle), nCycle = 1; end
if nargin < 4 || isempty(nOff), nOff = 0; end
if nargin < 5 || isempty(only), only = 1:nRuns; end
if nargin < 6 || isempty(duration), duration = 12600; end
if nargin < 7 || isempty(kScale), kScale = 1; end
if nargin < 8 || isempty(tauS), tauS = 0; end
p0 = smallsat_params;
orbit = 2*pi*sqrt((p0.Re + p0.altitude)^3/p0.mu);
mdl = 'smallsat_adcs_interference';
wasLoaded = bdIsLoaded(mdl);
load_system(mdl);
cleanup = onCleanup(@() evalin('base', 'clear adcsCase paramOverride'));
assignin('base', 'adcsCase', 'tumble');
rng(2026);
D = zeros(3, nRuns); R = zeros(1, nRuns); Q = zeros(4, nRuns);
for i = 1:nRuns                                   % same draws, same order as run_monte_carlo
    d = randn(3, 1); D(:, i) = d/norm(d);
    R(i) = (3 + 3*rand)*pi/180;
    qi = randn(4, 1); Q(:, i) = qi/norm(qi);
    rand(3, 1);
end
fprintf('Simulink copy: c = %.2e, cycle %d s, off %d s, tau %.2f s, gain x%.3f\n', c, nCycle, nOff, tauS, kScale);
res = struct('run', num2cell(1:nRuns));
for i = only
    ov = struct('tipOffRate', R(i)*D(:, i), 'initialQuat', Q(:, i), 'mtqLeak', c, 'mtqCycle', nCycle, ...
        'mtqOff', nOff, 'mtqLeakTau', tauS, 'bdotScale', kScale);
    assignin('base', 'paramOverride', ov);
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StopTime', num2str(duration));
    in = in.setBlockParameter([mdl '/Sensors/noise'], 'Seed', mat2str([11 22 33] + 1000*i));
    out = sim(in);
    t = out.t_log(:);
    tauW = log_rows(out.tauW_log);
    i0 = find(any(tauW ~= 0, 2), 1);          % hand-over to the wheels, as in run_monte_carlo
    if isempty(i0), tDet = NaN; else, tDet = t(i0); end
    wn = vecnorm(log_rows(out.w_log), 2, 2);
    m = log_rows(out.m_log);
    pre = t <= min(tDet, t(end));
    res(i).tip = R(i)*180/pi; res(i).tDet = tDet; res(i).wEnd = wn(end)*180/pi; res(i).mMax = max(max(abs(m(pre, :))));
    res(i).ok = ~isnan(tDet) && tDet <= 2*orbit;
    fprintf('%3d  tip-off %.2f deg/s  detumble %7.0f s  rate at end %.3f deg/s  peak dipole before hand-over %.1f  %s\n', ...
        i, res(i).tip, tDet, res(i).wEnd, res(i).mMax, ternary(res(i).ok, 'ok', 'MISS'));
end
res = res(only);
t = [res.tDet];
fprintf('%d of %d within 2 orbits (%.0f s); mean %.0f s over the %d runs that finished inside %.0f s (%d did not finish, so the mean is optimistic)\n', ...
    sum([res.ok]), numel(res), 2*orbit, mean(t(~isnan(t))), sum(~isnan(t)), duration, sum(isnan(t)));
if ~wasLoaded, close_system(mdl, 0); end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
