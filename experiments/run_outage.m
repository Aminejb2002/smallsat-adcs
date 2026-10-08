function res = run_outage(durations)
% Star tracker outage in the closed loop (estimate feeds the controller): the filter keeps
% propagating on the gyro and its torque model while no star fixes arrive. The outage starts at
% 14000 s, and lasts each of the given durations in seconds (0 = none).
% Needs the outage_on / outage_off Step blocks and the outage_valid gate in the Estimator; the gate is driven by
% P.outageStart and P.outageDuration (paramOverride).
% Reports the true knowledge error against the filter's own 3-sigma bound, and the pointing
% error during and after the outage. Logged every 10 s.
if nargin < 1, durations = [0 300 900 1800]; end
tStart = 14000;
p0 = smallsat_params;
arc = 180/pi*3600;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
cleanup = onCleanup(@() evalin('base', 'clear adcsCase paramOverride'));
res = struct('duration', num2cell(durations));
fprintf('%8s | %-34s | %-26s | %-24s\n', '', 'knowledge error in the outage', 'filter 3 sigma', 'pointing error [deg]');
fprintf('%8s | %10s %10s %10s | %12s %12s | %11s %11s %11s %9s\n', 'outage', 'rms', 'max', 'inside 3s', 'at end', 'before', 'max before', 'max during', 'max after', 'recovery');
fprintf('%8s | %10s %10s %10s | %12s %12s | %11s %11s %11s %9s\n', '[s]', '[arcsec]', '[arcsec]', '[%]', '[arcsec]', '[arcsec]', '-600 s', '', '+600 s', 'to <0.03');
for c = 1:numel(durations)
    dur = durations(c);
    if dur == 0, st = 1e9; else, st = tStart; end
    assignin('base', 'paramOverride', struct('outageStart', st, 'outageDuration', max(dur, 1)));
    in = Simulink.SimulationInput(mdl).setBlockParameter([mdl '/use_estimate'], 'Value', '1');
    out = sim(in);
    t = out.t_log(:);
    q = log_rows(out.q_log);
    qref = log_rows(out.qref_log);
    qHat = log_rows(out.estQ_log);
    sig = log_rows(out.estSig_log);
    idx = 1:round(10/p0.gyroStep):size(qHat, 1);
    qHat = qHat(idx, :);
    sig = sig(idx, :);
    n = numel(t);
    e = zeros(n, 3);
    pt = zeros(n, 1);
    for k = 1:n
        d = quat_mul([qHat(k, 1); -qHat(k, 2:4)'], q(k, :)');
        if d(1) < 0, d = -d; end
        e(k, :) = 2*d(2:4)';
        qe = quat_mul([qref(k, 1); -qref(k, 2:4)'], q(k, :)');
        pt(k) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
    end
    win = t >= tStart & t < tStart + max(dur, 1000);
    before = t >= tStart - 600 & t < tStart;
    after = t >= tStart + dur & t < tStart + dur + 600;
    inside = mean(all(abs(e(win, :)) < 3*sig(win, :), 2))*100;
    eN = vecnorm(e, 2, 2)*arc;
    sN = vecnorm(sig, 2, 2)*3*arc;
    iEnd = find(t >= tStart + dur, 1) - 1;   % last logged sample inside the outage, before the first star update
    res(c).rms = sqrt(mean(sum(e(win, :).^2, 2)))*arc;   % total angle, same definition as run_monte_carlo
    res(c).max = max(eN(win));
    res(c).inside = inside;
    res(c).sig3End = sN(iEnd);
    res(c).sig3Before = mean(sN(before));
    late = find(t >= tStart + dur & pt >= 0.03, 1, 'last');
    if isempty(late), res(c).recover = 0; else, res(c).recover = t(late) + 10 - (tStart + dur); end
    res(c).pointMaxBefore = max(pt(before));
    res(c).pointMaxDuring = max(pt(win));
    res(c).pointMaxAfter = max(pt(after));
    fprintf('%8.0f | %10.1f %10.1f %10.1f | %12.1f %12.1f | %11.4f %11.4f %11.4f %9.0f\n', dur, res(c).rms, res(c).max, ...
        inside, res(c).sig3End, res(c).sig3Before, res(c).pointMaxBefore, res(c).pointMaxDuring, res(c).pointMaxAfter, res(c).recover);
end
if ~wasLoaded, close_system(mdl, 0); end
end
