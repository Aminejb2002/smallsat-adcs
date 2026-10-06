function estimator_crosscheck
% Runs models/smallsat_adcs.slx for 600 s with the sensor noise on, then feeds the logged filter
% inputs (rate, star quaternion, update flag) through src/estimation again in plain MATLAB and
% compares the estimates. Two cases: the 'point' start, where star fixes arrive and the filter
% runs, and the 'tumble' start, where the rate is above P.starRateLimit and the filter must
% stay silent and pass the raw measurements through. Needs these logs at P.gyroStep:
% estIn_w_log, estIn_q_log, estIn_new_log, estQ_log, estB_log, estSig_log.
% Thresholds were fixed before the first run.
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
cleanup = onCleanup(@() evalin('base', 'clear adcsCase'));
duration = 600;
expectedSamples = round(duration/p.gyroStep) + 1;
cases = {'point', floor(duration/p.starStep) + 1; 'tumble', 0};
fails = 0;
total = 0;
for c = 1:size(cases, 1)
    assignin('base', 'adcsCase', cases{c, 1});
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StopTime', num2str(duration));
    out = sim(in);

    w = log_rows(out.estIn_w_log);
    qm = log_rows(out.estIn_q_log);
    new = log_rows(out.estIn_new_log);
    qSim = log_rows(out.estQ_log);
    bSim = log_rows(out.estB_log);
    sigSim = log_rows(out.estSig_log);
    n = size(w, 1);

    x = zeros(43, 1);
    qRef = zeros(n, 4);
    bRef = zeros(n, 3);
    sigRef = zeros(n, 3);
    for k = 1:n
        x = mekf_step(x, w(k, :)', qm(k, :)', new(k), p);
        [qh, ~, sg, bh] = mekf_outputs(x, w(k, :)', qm(k, :)');
        qRef(k, :) = qh';
        bRef(k, :) = bh';
        sigRef(k, :) = sg';
    end

    sgn = sign(sum(qSim.*qRef, 2));
    sgn(sgn == 0) = 1;
    checks = {
        'sample count',        abs(n - expectedSamples),                          0.5
        'update count',        abs(nnz(new) - cases{c, 2}),                       0.5
        'attitude estimate',   max(2*vecnorm(qSim.*sgn - qRef, 2, 2)),            1e-9
        'bias estimate',       max(vecnorm(bSim - bRef, 2, 2)),                   1e-12
        'attitude sigma',      max(vecnorm(sigSim - sigRef, 2, 2)),               1e-12
        };
    fprintf('\n%s start: %d samples, %d star updates (expected %d)\n', cases{c, 1}, n, nnz(new), cases{c, 2});
    for k = 1:size(checks, 1)
        ok = checks{k, 2} < checks{k, 3};
        if ok, s = 'PASS'; else, s = 'FAIL'; fails = fails + 1; end
        total = total + 1;
        fprintf('%s  %-20s max difference %.3e (limit %.0e)\n', s, checks{k, 1}, checks{k, 2}, checks{k, 3});
    end
end
if ~wasLoaded
    close_system(mdl, 0);
end
fprintf('\n%d of %d checks passed\n', total - fails, total);
if fails > 0
    error('estimator_crosscheck:failed', '%d checks failed', fails);
end
end
