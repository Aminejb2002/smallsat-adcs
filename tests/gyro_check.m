function gyro_check
% Statistical check of the gyro model in models/smallsat_adcs.slx over 6000 s. Needs the logs
% gyroErr_log (measured rate minus true rate) and gyroBias_log (bias state), sampled every
% P.gyroStep. Thresholds were fixed before the first run.
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
duration = 6000;
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StopTime', num2str(duration));
out = sim(in);
if ~wasLoaded
    close_system(mdl, 0);
end

err = log_rows(out.gyroErr_log);
bias = log_rows(out.gyroBias_log);
noise = err - bias;
n = size(err, 1);
sigma = p.gyroArw/sqrt(p.gyroStep);
driftSigma = p.gyroRrw*sqrt(duration);

checks = {
    'initial bias equals gyroBias0',   max(abs(bias(1, :)' - p.gyroBias0)),                   1e-12
    'noise mean',                      max(abs(mean(noise, 1)))/(sigma/sqrt(n)),                4
    'noise std vs ARW/sqrt(Ts)',       max(abs(std(noise, 0, 1)/sigma - 1)),                    0.03
    'bias step std vs RRW*sqrt(Ts)',   max(abs(std(diff(bias, 1, 1), 0, 1)/(p.gyroRrw*sqrt(p.gyroStep)) - 1)), 0.03
    'bias end drift (not frozen)',     max(abs(bias(end, :) - bias(1, :)))/driftSigma,          4
    };
fails = 0;
for k = 1:size(checks, 1)
    ok = checks{k, 2} < checks{k, 3};
    if ok, s = 'PASS'; else, s = 'FAIL'; fails = fails + 1; end
    fprintf('%s  %-32s %10.3e (limit %.0e)\n', s, checks{k, 1}, checks{k, 2}, checks{k, 3});
end
fprintf('%d samples, noise std %.3e rad/s (expected %.3e)\n', n, mean(std(noise, 0, 1)), sigma);
fprintf('bias step std %.3e rad/s (expected %.3e), end drift %.3e rad/s (expected scale %.3e)\n', ...
    mean(std(diff(bias, 1, 1), 0, 1)), p.gyroRrw*sqrt(p.gyroStep), mean(abs(bias(end, :) - bias(1, :))), driftSigma);
fprintf('\n%d of %d checks passed\n', size(checks, 1) - fails, size(checks, 1));
if fails > 0
    error('gyro_check:failed', '%d checks failed', fails);
end
end
