function estimator_crosscheck
% Runs models/smallsat_adcs.slx for 600 s with the sensor noise on, then feeds the logged filter
% inputs (rate, star quaternion, update pulse) through src/estimation again in plain MATLAB and
% compares the estimates. Needs these logs at P.gyroStep: estIn_w_log, estIn_q_log,
% estIn_new_log, estQ_log, estB_log, estSig_log. Thresholds were fixed before the first run.
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
duration = 600;
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StopTime', num2str(duration));
out = sim(in);
if ~wasLoaded
    close_system(mdl, 0);
end

w = log_rows(out.estIn_w_log);
qm = log_rows(out.estIn_q_log);
new = log_rows(out.estIn_new_log);
qSim = log_rows(out.estQ_log);
bSim = log_rows(out.estB_log);
sigSim = log_rows(out.estSig_log);
n = size(w, 1);
expected = round(duration/p.gyroStep) + 1;
fprintf('%d samples (expected %d), %d star updates\n', n, expected, nnz(new));

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
checks = {
    'sample count',            abs(n - expected),                                   0.5
    'update pulse count',      abs(nnz(new) - (floor(duration/p.starStep) + 1)),   0.5
    'attitude estimate',       max(2*vecnorm(qSim.*sgn - qRef, 2, 2)),              1e-9
    'bias estimate',           max(vecnorm(bSim - bRef, 2, 2)),                     1e-12
    'attitude sigma',          max(vecnorm(sigSim - sigRef, 2, 2)),                 1e-12
    };
fails = 0;
for k = 1:size(checks, 1)
    ok = checks{k, 2} < checks{k, 3};
    if ok, s = 'PASS'; else, s = 'FAIL'; fails = fails + 1; end
    fprintf('%s  %-22s max difference %.3e (limit %.0e)\n', s, checks{k, 1}, checks{k, 2}, checks{k, 3});
end
fprintf('\n%d of %d checks passed\n', size(checks, 1) - fails, size(checks, 1));
if fails > 0
    error('estimator_crosscheck:failed', '%d checks failed', fails);
end
end
