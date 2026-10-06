function res = run_pointing_sweep
% Pointing loop bandwidth and damping against the gyro noise that reaches the wheel torque
% through the rate term. Three-orbit tumble case with the estimate in the loop, one run per
% (wn, zeta). Each run uses the same noise seeds, so differences between rows come from the
% gains only; the tail statistics come from one noise realisation. Logged quantities are
% sampled every 10 s, so torque peaks between samples are not seen.
grid = [0.05 0.9; 0.05 0.7; 0.075 0.7; 0.10 0.9; 0.10 0.7; 0.15 0.7; 0.20 0.7];
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
cleanup = onCleanup(@() evalin('base', 'clear adcsCase paramOverride'));
nCase = size(grid, 1);
res = struct('wn', num2cell(grid(:, 1)), 'zeta', num2cell(grid(:, 2)));
fprintf('%6s %5s %9s %9s %9s %11s %11s %9s %9s\n', 'wn', 'zeta', 'kd', 'handover', 'capture', 'mean err', 'max err', 'sat', 'h peak');
fprintf('%6s %5s %9s %9s %9s %11s %11s %9s %9s\n', '[rad/s]', '', '[Nms/rad]', '[s]', '[s]', '[deg]', '[deg]', 'fraction', '[Nms]');
for c = 1:nCase
    assignin('base', 'paramOverride', struct('pointWn', grid(c, 1), 'pointZeta', grid(c, 2)));
    p = smallsat_params;
    in = Simulink.SimulationInput(mdl).setBlockParameter([mdl '/use_estimate'], 'Value', '1');
    out = sim(in);
    t = out.t_log(:);
    q = log_rows(out.q_log);
    qref = log_rows(out.qref_log);
    tauW = log_rows(out.tauW_log);
    h = log_rows(out.h_log);
    errDeg = zeros(numel(t), 1);
    for k = 1:numel(t)
        qe = quat_mul([qref(k, 1); -qref(k, 2:4)'], q(k, :)');
        errDeg(k) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
    end
    i0 = find(any(tauW ~= 0, 2), 1);
    tEngage = t(i0);
    bad = find(errDeg >= 0.1 & t >= tEngage, 1, 'last');
    if isempty(bad)
        capture = 0;
    elseif bad == numel(t)
        capture = NaN;
    else
        capture = t(bad + 1) - tEngage;
    end
    tail = t >= t(end) - 2000;
    sat = mean(any(abs(tauW(t >= tEngage, :)) >= p.wheelTorqueMax - 1e-9, 2));
    res(c).handover = tEngage;
    res(c).capture = capture;
    res(c).errMean = mean(errDeg(tail));
    res(c).errMax = max(errDeg(tail));
    res(c).sat = sat;
    res(c).hPeak = max(vecnorm(h, 2, 2));
    kd = 2*p.pointZeta*p.pointWn*p.inertia(3, 3);
    fprintf('%6.3f %5.2f %9.2f %9.0f %9.0f %11.4f %11.4f %9.3f %9.3f\n', grid(c, 1), grid(c, 2), kd, ...
        tEngage, capture, res(c).errMean, res(c).errMax, sat, res(c).hPeak);
end
if ~wasLoaded
    close_system(mdl, 0);
end
end
