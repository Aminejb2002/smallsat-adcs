function res = run_estimator_robustness
% Accuracy of the rate-aided filter when what it is told differs from the truth: inertia error,
% wheel torque error, an unmodelled wheel drag torque, torque noise and a star tracker that is
% worse about its boresight (body z); the filter knows ratio 5 unless a case says otherwise. Filter only, no controller:
% the attitude error here is the knowledge error, which sets the floor of the pointing error.
% Errors are rms per axis over t >= 600 s, mean over seeds. NEES uses the filter's own covariance.
% A weak rate feedback (rateDamp, as a controller would) keeps the body rate bounded; the 'limit' row
% switches it off, so the drag spins the body up to about 1 deg/s and the filter loses consistency.
% The wheel drag, torque noise and boresight noise values are assumptions, not data of a real unit.
p = smallsat_params;
arcsec = 180/pi*3600;
sigGyro = p.gyroArw/sqrt(p.gyroStep)*180/pi;
cases = {
    'ideal knowledge',            struct()
    'inertia +10 %',              struct('inertiaErr', 0.10)
    'inertia +10/-10/+5 %',       struct('inertiaErr', [0.10; -0.10; 0.05])
    'wheel torque gain +5 %',     struct('torqueErr', 0.05)
    'wheel drag 2e-4 N m',        struct('torqueBias', [2e-4; -2e-4; 2e-4])
    'torque noise 2e-4 N m',      struct('torqueNoise', 2e-4)
    'realistic nominal',          struct('inertiaErr', [0.10; -0.10; 0.05], 'torqueErr', 0.05, ...
                                         'torqueBias', [2e-4; -2e-4; 2e-4])
    'limit: drag, no controller', struct('inertiaErr', [0.10; -0.10; 0.05], 'torqueBias', [2e-4; -2e-4; 2e-4], 'rateDamp', 0)
    'boresight ratio 5, unaware', struct('starRatioFilter', 1)
    'boresight ratio 10, filter 5', struct('starRatioTrue', 10)
    'stress: all of the above',           struct('inertiaErr', [0.10; -0.10; 0.05], 'torqueErr', 0.05, ...
                                         'torqueBias', [2e-4; -2e-4; 2e-4], 'torqueNoise', 2e-4, ...
                                         'starRatioTrue', 10)
    };
seeds = [11 12];
T = 3000;
nCase = size(cases, 1);
res = struct('name', cases(:, 1), 'att', 0, 'attMax', 0, 'rate', 0, 'neesAtt', 0, 'neesRate', 0, 'bias', 0);
fprintf('%-26s %10s %10s %12s %12s %9s %9s\n', 'case', 'att rms', 'att max', 'rate rms', 'rate / raw', 'NEES att', 'NEES rate');
fprintf('%-26s %10s %10s %12s %12s %9s %9s\n', '', '[arcsec]', '[arcsec]', '[deg/s]', 'gyro noise', '(ideal 3)', '(ideal 3)');
for c = 1:nCase
    att = 0; attMax = 0; rate = 0; na = 0; nr = 0;
    for s = 1:numel(seeds)
        rec = estimator_scenario(p, seeds(s), T, [], true, p.gyroBias0, 0, cases{c, 2});
        sel = find(rec.t >= 600);
        e = rec.e(:, sel);
        wE = rec.wErr(:, sel);
        att = att + sqrt(mean(e(:).^2))*arcsec/numel(seeds);
        attMax = max(attMax, max(vecnorm(e, 2, 1))*arcsec);
        rate = rate + sqrt(mean(wE(:).^2))*180/pi/numel(seeds);
        a = zeros(numel(sel), 1);
        r = zeros(numel(sel), 1);
        for k = 1:numel(sel)
            a(k) = rec.e(:, sel(k))'*(rec.P3(:, :, sel(k)) \ rec.e(:, sel(k)));
            r(k) = rec.wErr(:, sel(k))'*(rec.Pw(:, :, sel(k)) \ rec.wErr(:, sel(k)));
        end
        na = na + mean(a)/numel(seeds);
        nr = nr + mean(r)/numel(seeds);
    end
    res(c).att = att;
    res(c).attMax = attMax;
    res(c).rate = rate;
    res(c).neesAtt = na;
    res(c).neesRate = nr;
    fprintf('%-26s %10.1f %10.1f %12.2e %12.3f %9.2f %9.2f\n', cases{c, 1}, att, attMax, rate, rate/sigGyro, na, nr);
end
fprintf('\nreference, raw gyro rate minus bias: %.2e deg/s per axis; star tracker %.1f arcsec per axis\n', sigGyro, p.starSigma*arcsec);
end
