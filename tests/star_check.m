function star_check
% Statistical check of the star tracker model in models/smallsat_adcs.slx over 6000 s. Needs the
% logs starTrue_log (held true quaternion) and starMeas_log (measured quaternion), sampled every
% P.starStep. The error is the small rotation e with qMeas = qTrue (x) [1; e/2], taken in the
% body frame. Thresholds were fixed before the first run.
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

qt = log_rows(out.starTrue_log);
qm = log_rows(out.starMeas_log);
n = size(qt, 1);
e = zeros(n, 3);
for k = 1:n
    d = quat_mul([qt(k, 1); -qt(k, 2:4)'], qm(k, :)');
    if d(1) < 0, d = -d; end
    e(k, :) = 2*d(2:4)';
end
sigmaAx = sqrt(diag(star_covariance(p)))';
sigma = p.starSigma;

checks = {
    'measured quaternion norm',        max(abs(vecnorm(qm, 2, 2) - 1)),                1e-9
    'error mean',                      max(abs(mean(e, 1))./(sigmaAx/sqrt(n))),            4
    'error std vs star covariance',    max(abs(std(e, 0, 1)./sigmaAx - 1)),                0.05
    'axes uncorrelated',               max(abs(corrcoef_offdiag(e))),                   4/sqrt(n)
    };
fails = 0;
for k = 1:size(checks, 1)
    ok = checks{k, 2} < checks{k, 3};
    if ok, s = 'PASS'; else, s = 'FAIL'; fails = fails + 1; end
    fprintf('%s  %-28s %10.3e (limit %.1e)\n', s, checks{k, 1}, checks{k, 2}, checks{k, 3});
end
fprintf('%d samples, error std %.3e rad per axis (expected %.3e, %.1f arcsec)\n', n, ...
    mean(std(e, 0, 1)), sigma, sigma*180/pi*3600);
fprintf('\n%d of %d checks passed\n', size(checks, 1) - fails, size(checks, 1));
if fails > 0
    error('star_check:failed', '%d checks failed', fails);
end
end

function c = corrcoef_offdiag(e)
r = corrcoef(e);
c = [r(1, 2) r(1, 3) r(2, 3)];
end
