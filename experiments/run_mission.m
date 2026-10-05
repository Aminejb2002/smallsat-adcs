function out = run_mission
% First end-to-end run: tip-off tumble, B-dot detumble, handover, wheel capture of nadir
% pointing, three orbits, noisy magnetometer. Truth state feeds the controller (no estimator yet).
% Quantities are sampled every 10 s, so peaks between samples are not seen.
root = fileparts(fileparts(mfilename('fullpath')));
p = smallsat_params;
bd = size_bdot(p);
pt = size_pointing(p);
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
simOut = sim(mdl);
evalin('base', 'clear adcsCase');
if ~wasLoaded
    close_system(mdl, 0);
end

t = simOut.t_log(:);
w = log_rows(simOut.w_log);
q = log_rows(simOut.q_log);
qref = log_rows(simOut.qref_log);
h = log_rows(simOut.h_log);
m = log_rows(simOut.m_log);
tauW = log_rows(simOut.tauW_log);
period = 2*pi*sqrt((p.Re + p.altitude)^3/p.mu);

rate = vecnorm(w, 2, 2);
errDeg = zeros(numel(t), 1);
for k = 1:numel(t)
    qe = quat_mul([qref(k, 1); -qref(k, 2:4)'], q(k, :)');
    errDeg(k) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
end

engaged = find(any(tauW ~= 0, 2), 1);
if isempty(engaged)
    error('run_mission:noHandover', 'wheels never engaged within %.0f s', t(end));
end
tEngage = t(engaged);
bad = find(errDeg >= 0.1 & t >= tEngage, 1, 'last');
if isempty(bad)
    tCapture = tEngage;
elseif bad == numel(t)
    tCapture = NaN;
else
    tCapture = t(bad + 1);
end
tail = t >= t(end) - 2000;
torquersOff = all(all(m(t >= tEngage + 10, :) == 0));
engagedRows = t >= tEngage;
satFraction = mean(any(abs(tauW(engagedRows, :)) >= p.wheelTorqueMax - 1e-9, 2));

out.t = t;
out.rate = rate;
out.errDeg = errDeg;
out.tEngage = tEngage;
out.tCapture = tCapture;
out.hPeak = max(vecnorm(h, 2, 2));

fprintf('tip-off %.2f deg/s, handover below %.2f deg/s, exit above %.2f deg/s\n', ...
    norm(p.tipOffRate)*180/pi, p.detumbleTarget*180/pi, p.detumbleExit*180/pi);
fprintf('%-38s %12s %12s\n', '', 'hand sizing', 'simulation');
fprintf('%-38s %12.0f %12.0f\n', 'B-dot time to handover rate [s]', bd.tDetumble, tEngage);
fprintf('%-38s %12.0f %12.0f\n', 'wheel capture, start to 0.1 deg [s]', pt.tCapture, tCapture - tEngage);
fprintf('%-38s %12.2f %12.2f\n', 'peak wheel momentum [N m s]', pt.hSlew, out.hPeak);
fprintf('(hand capture time is the 180 deg worst case; hand momentum is the slew-rate-limit value)\n');
fprintf('torquers silent after handover: %s\n', yesno(torquersOff));
fprintf('wheel torque at its limit in %.1f %% of samples after handover\n', 100*satFraction);
fprintf('pointing error over the last %.0f s: mean %.4f deg, max %.4f deg\n', 2000, ...
    mean(errDeg(tail)), max(errDeg(tail)));
fprintf('final rate %.4f deg/s, final wheel momentum %.3f N m s (limit %.1f)\n', ...
    rate(end)*180/pi, norm(h(end, :)), p.wheelMomentumMax);

fig = figure('Name', 'Mission run', 'Position', [100 100 900 800]);
subplot(3, 1, 1);
semilogy(t/period, rate*180/pi, 'LineWidth', 1.2);
hold on;
yline(p.detumbleTarget*180/pi, '--');
yline(p.detumbleExit*180/pi, ':');
grid on;
ylabel('|\omega| [deg/s]');
title('Detumble, handover and pointing');
subplot(3, 1, 2);
semilogy(t/period, max(errDeg, 1e-4), 'LineWidth', 1.2);
hold on;
yline(0.1, '--');
grid on;
ylabel('pointing error [deg]');
subplot(3, 1, 3);
plot(t/period, h, 'LineWidth', 1.0);
hold on;
plot(t/period, vecnorm(h, 2, 2), 'k', 'LineWidth', 1.2);
grid on;
xlabel('orbits');
ylabel('wheel momentum [N m s]');
legend('x', 'y', 'z', 'norm', 'Location', 'best');
mediaDir = fullfile(root, 'results', 'media');
if ~exist(mediaDir, 'dir'), mkdir(mediaDir); end
exportgraphics(fig, fullfile(mediaDir, 'mission_run_draft.png'), 'Resolution', 150);

if nargout == 0
    clear out
end
end

function s = yesno(flag)
if flag
    s = 'yes';
else
    s = 'NO';
end
end
