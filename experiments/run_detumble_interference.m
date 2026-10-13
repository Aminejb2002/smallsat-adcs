function res = run_detumble_interference(nRuns, c, nCycle, nOff, only, duration, kScale, tauS)
% B-dot detumble with magnetorquer field leaking into the magnetometer, plain MATLAB.
% This is NOT the Simulink model: it uses the same sampled B-dot loop as reference_run (orbit,
% attitude, disturbances, tumble cases) with a noisy magnetometer, so the baseline (c = 0) can be
% compared with the Simulink Monte Carlo before the interference is switched on.
%
%   c       coupling, T per A m^2: B_meas = B + c*m_applied, same axis, instantaneous (assumption).
%           Dipole-on-axis estimate: c ~ 2e-7/d^3 for coil-to-sensor distance d in m.
%   nCycle  samples per B-dot cycle in seconds (1 = continuous, torquers always on while sampling)
%   nOff    seconds with torquers off at the start of each cycle. The field is read at the end of the
%           off window (settled, no leakage), the new dipole is held for the other nCycle - nOff seconds.
%           Duty cycle = (nCycle - nOff)/nCycle.
%   only    run numbers to simulate (tip-off, attitude and gyro-bias draws are the same as in
%           run_monte_carlo, so run i has the same initial condition; the magnetometer noise differs)
%   duration  s, default 2 orbits (11 600 s) plus 1000 s so a late finish is still seen
%   tauS    settling time constant of the leakage after the torquers switch off, s (default 0 = instant).
%           At the read, the leakage left is c*m_last*exp(-nOff/tauS). With nOff = 0 it is the full c*m.
%   kScale  factor on the B-dot gain from size_bdot (default 1); e.g. 1/duty to make up for the off-time
% Detumble time: first time the rate is below p.detumbleTarget and stays below it.
if nargin < 1 || isempty(nRuns), nRuns = 30; end
if nargin < 2 || isempty(c), c = 0; end
if nargin < 3 || isempty(nCycle), nCycle = 1; end
if nargin < 4 || isempty(nOff), nOff = 0; end
if nargin < 5 || isempty(only), only = 1:nRuns; end
if nargin < 6 || isempty(duration), duration = 12600; end
if nargin < 7 || isempty(kScale), kScale = 1; end
if nargin < 8 || isempty(tauS), tauS = 0; end
p = smallsat_params;
bd = size_bdot(p);
bd.k = bd.k*kScale;
Ts = p.sensorStep;
nStep = round(duration/Ts);
orbit = 2*pi*sqrt((p.Re + p.altitude)^3/p.mu);
rng(2026);
res = struct('run', num2cell(1:nRuns));
fprintf('c = %.2e T/(A m^2), cycle %d s, off %d s (duty %.0f %%), gain k = %.3e, loop gain k*c/cycle = %.2f\n', ...
    c, nCycle, nOff, 100*(nCycle - nOff)/nCycle, bd.k, bd.k*c/nCycle);
D = zeros(3, nRuns); R = zeros(1, nRuns); Q = zeros(4, nRuns);
for i = 1:nRuns                                   % same draws, same order as run_monte_carlo
    d = randn(3, 1); D(:, i) = d/norm(d);
    R(i) = (3 + 3*rand)*pi/180;
    qi = randn(4, 1); Q(:, i) = qi/norm(qi);
    rand(3, 1);                                   % gyro bias, unused here
end
for i = only
    ic = initial_state_tumble(p);
    rate = R(i);
    ic.q = Q(:, i); ic.w = rate*D(:, i); ic.h = zeros(3, 1);
    rng(5000 + i);                                % magnetometer noise stream of this run
    x = [ic.r; ic.v; ic.q; ic.w; ic.h];
    wn = zeros(nStep + 1, 1); wn(1) = norm(x(11:13));
    mMax = 0; mHold = zeros(3, 1); mApplied = zeros(3, 1);
    Bold = [];
    for k = 0:nStep-1
        tk = k*Ts;
        phase = mod(k, nCycle);
        if phase == mod(nOff, nCycle)      % end of the off window: read the field
            Btrue = body_field(tk, x, p);
            if nOff > 0 && tauS > 0, leak = c*mHold*exp(-nOff/tauS);   % torquers were off for nOff s
            elseif nOff > 0, leak = zeros(3, 1);
            else, leak = c*mApplied; end
            Bm = Btrue + leak + p.magSigma*randn(3, 1);
            if isempty(Bold), Bold = Bm; end
            mHold = -bd.k*(Bm - Bold)/nCycle;
            mHold = min(max(mHold, -p.mtqMax), p.mtqMax);
            Bold = Bm;
        end
        if phase < nOff, mApplied = zeros(3, 1); else, mApplied = mHold; end
        mMax = max(mMax, max(abs(mApplied)));
        [~, y] = ode113(@(tt, xx) rhs(tt, xx, p, mApplied), [tk, tk + Ts], x, odeset('RelTol', 1e-9, 'AbsTol', 1e-9));
        x = y(end, :)';
        wn(k + 2) = norm(x(11:13));
    end
    last = find(wn >= p.detumbleTarget, 1, 'last');
    if isempty(last), tDet = 0; elseif last == numel(wn), tDet = NaN; else, tDet = last*Ts; end
    res(i).tip = rate*180/pi; res(i).tDet = tDet; res(i).mMax = mMax; res(i).wEnd = wn(end)*180/pi;
    res(i).ok = tDet <= 2*orbit;
    fprintf('%3d  tip-off %.2f deg/s  detumble %7.0f s  final %.3f deg/s  peak dipole %.1f  %s\n', ...
        i, res(i).tip, tDet, res(i).wEnd, mMax, ternary(res(i).ok, 'ok', 'MISS'));
end
res = res(only);
t = [res.tDet];
fprintf('%d of %d within 2 orbits (%.0f s); mean %.0f s over the %d runs that finished inside %.0f s (%d did not finish, so the mean is optimistic)\n', ...
    sum([res.ok]), numel(res), 2*orbit, mean(t(~isnan(t))), sum(~isnan(t)), duration, sum(isnan(t)));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function B = body_field(t, x, p)
B = quat_to_dcm(x(7:10))*magnetic_field_eci(x(1:3), t, p);
end

function dx = rhs(t, x, p, m)
tauDist = disturbance_block(t, x(1:3), x(4:6), x(7:10), p);
tauMt = cross(m, body_field(t, x, p));
wd = angular_accel(x(11:13), x(14:16), tauDist + tauMt, zeros(3, 1), p.inertiaTruth);
dx = [x(4:6); j2_accel(x(1:3), p); quat_rate(x(7:10), x(11:13)); wd; zeros(3, 1)];
end
