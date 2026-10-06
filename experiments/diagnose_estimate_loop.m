function diagnose_estimate_loop
% Three-orbit runs that separate what the estimate in the loop changes. Cases:
%   truth       use_estimate = 0, the truth-fed reference run
%   estimate    use_estimate = 1, everything on
%   no gyro noise  use_estimate = 1 with the gyro white noise gain set to zero
% For each case: handover time, true body rate at handover, wheel capture time and the
% pointing error over the last 2000 s (true attitude against the reference).
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
cleanup = onCleanup(@() evalin('base', 'clear adcsCase'));
base = Simulink.SimulationInput(mdl);

cases = {'truth', 0, []; 'estimate', 1, []; 'no gyro noise', 1, '0'};
fprintf('%-16s %10s %14s %12s %14s %14s\n', 'case', 'handover', 'true rate', 'capture', 'last 2000 s', 'last 2000 s');
fprintf('%-16s %10s %14s %12s %14s %14s\n', '', '[s]', '[deg/s]', '[s]', 'mean [deg]', 'max [deg]');
for c = 1:size(cases, 1)
    in = base.setBlockParameter([mdl '/use_estimate'], 'Value', num2str(cases{c, 2}));
    if ~isempty(cases{c, 3})
        in = in.setBlockParameter([mdl '/Sensors/noise_gain_w'], 'Gain', cases{c, 3});
    end
    out = sim(in);
    m = metrics(out, p);
    fprintf('%-16s %10.0f %14.3f %12.0f %14.4f %14.4f\n', cases{c, 1}, m.tEngage, m.rateAtEngage, ...
        m.capture, m.errMean, m.errMax);
end
if ~wasLoaded
    close_system(mdl, 0);
end
end

function m = metrics(out, p)
t = out.t_log(:);
w = log_rows(out.w_log);
q = log_rows(out.q_log);
qref = log_rows(out.qref_log);
tauW = log_rows(out.tauW_log);
errDeg = zeros(numel(t), 1);
for k = 1:numel(t)
    qe = quat_mul([qref(k, 1); -qref(k, 2:4)'], q(k, :)');
    errDeg(k) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
end
i0 = find(any(tauW ~= 0, 2), 1);
m.tEngage = t(i0);
m.rateAtEngage = norm(w(i0, :))*180/pi;
bad = find(errDeg >= 0.1 & t >= m.tEngage, 1, 'last');
if isempty(bad)
    m.capture = 0;
elseif bad == numel(t)
    m.capture = NaN;
else
    m.capture = t(bad + 1) - m.tEngage;
end
tail = t >= t(end) - 2000;
m.errMean = mean(errDeg(tail));
m.errMax = max(errDeg(tail));
end
