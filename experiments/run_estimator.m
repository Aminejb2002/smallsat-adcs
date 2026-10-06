function out = run_estimator
% Full three-orbit run (tumble, detumble, pointing) with the MEKF running on the noisy gyro and
% star tracker. The controller still uses the true state; this run only judges the estimate.
% Error is the body-frame rotation between estimate and truth. Needs the Estimator logs and
% gyroBias_log in the model. Truth is logged every 10 s, estimates every P.gyroStep.
root = fileparts(fileparts(mfilename('fullpath')));
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
simOut = sim(mdl);
evalin('base', 'clear adcsCase');
if ~wasLoaded
    close_system(mdl, 0);
end

step = round(10/p.gyroStep);
t = simOut.t_log(:);
qTrue = log_rows(simOut.q_log);
w = log_rows(simOut.w_log);
qHat = log_rows(simOut.estQ_log);
sig = log_rows(simOut.estSig_log);
bHat = log_rows(simOut.estB_log);
bTrue = log_rows(simOut.gyroBias_log);
idx = 1:step:size(qHat, 1);
if numel(idx) ~= numel(t)
    error('run_estimator:grid', 'estimate log has %d samples at 10 s, truth has %d', numel(idx), numel(t));
end
qHat = qHat(idx, :);
sig = sig(idx, :);
bHat = bHat(idx, :);
bTrue = bTrue(idx, :);

n = numel(t);
e = zeros(n, 3);
for k = 1:n
    d = quat_mul([qHat(k, 1); -qHat(k, 2:4)'], qTrue(k, :)');
    if d(1) < 0, d = -d; end
    e(k, :) = 2*d(2:4)';
end
period = 2*pi*sqrt((p.Re + p.altitude)^3/p.mu);
active = sig(:, 1) > 0;
tInit = t(find(active, 1));
if isempty(tInit)
    error('run_estimator:noFix', 'the filter never started: no valid star tracker fix in %.0f s', t(end));
end
iInit = find(active, 1);
settled = active & t >= tInit + 600;
last = active & t >= t(end) - period;
e(~active, :) = NaN;
arcsec = 180/pi*3600;
rmsAll = sqrt(mean(sum(e(settled, :).^2, 2)/3));
rmsLast = sqrt(mean(sum(e(last, :).^2, 2)/3));
inside = mean(all(abs(e(settled, :)) < 3*sig(settled, :), 2));
bErr = vecnorm(bTrue - bHat, 2, 2);
bErr(~active) = NaN;

out.t = t;
out.err = e;
out.sig = sig;
out.bErr = bErr;

fprintf('filter starts at %.0f s (first star fix, body rate %.2f deg/s, limit %.1f deg/s), statistics from %.0f s\n', ...
    tInit, norm(w(iInit, :))*180/pi, p.starRateLimit*180/pi, tInit + 600);
fprintf('attitude error per axis, rms after settling: %.2f arcsec (star tracker alone %.1f arcsec)\n', ...
    rmsAll*arcsec, p.starSigma*arcsec);
fprintf('attitude error per axis, rms in the last orbit: %.2f arcsec\n', rmsLast*arcsec);
fprintf('samples inside the filter 3 sigma bound: %.1f %%\n', 100*inside);
fprintf('largest error after settling: %.1f arcsec (during rates up to %.2f deg/s)\n', ...
    max(vecnorm(e(settled, :), 2, 2))*arcsec, max(vecnorm(w(settled, :), 2, 2))*180/pi);
fprintf('bias error: at start of filter %.2e rad/s, end %.2e rad/s, true bias norm %.2e rad/s\n', ...
    bErr(iInit), bErr(end), norm(bTrue(end, :)));

fig = figure('Name', 'Estimator run', 'Position', [100 100 900 800]);
subplot(3, 1, 1);
sigPlot = sig;
sigPlot(~active, :) = NaN;
semilogy(t/period, vecnorm(e, 2, 2)*arcsec, 'LineWidth', 1.0);
hold on;
semilogy(t/period, vecnorm(sigPlot, 2, 2)*arcsec*3, '--', 'LineWidth', 1.0);
grid on;
ylabel('attitude error [arcsec]');
legend('error norm', '3 sigma bound (norm)', 'Location', 'best');
title('MEKF on gyro and star tracker');
subplot(3, 1, 2);
plot(t/period, (bTrue - bHat)*180/pi*3600, 'LineWidth', 1.0);
grid on;
ylabel('bias error [deg/h]');
subplot(3, 1, 3);
plot(t/period, vecnorm(w, 2, 2)*180/pi, 'LineWidth', 1.0);
grid on;
xlabel('orbits');
ylabel('|\omega| [deg/s]');
mediaDir = fullfile(root, 'results', 'media');
if ~exist(mediaDir, 'dir'), mkdir(mediaDir); end
exportgraphics(fig, fullfile(mediaDir, 'estimator_run_draft.png'), 'Resolution', 150);
if nargout == 0
    clear out
end
end
