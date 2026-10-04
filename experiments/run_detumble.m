function out = run_detumble
% B-dot detumble from the tip-off rate, three orbits, noisy magnetometer. Compares the
% simulation with the hand sizing in size_bdot.
p = smallsat_params;
bd = size_bdot(p);
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
m = log_rows(simOut.m_log);
wn = vecnorm(w, 2, 2);
period = 2*pi*sqrt((p.Re + p.altitude)^3/p.mu);

last = find(wn >= p.detumbleTarget, 1, 'last');
if isempty(last)
    tTarget = 0;
elseif last == numel(t)
    tTarget = NaN;
else
    tTarget = t(last + 1);
end

window = wn < 0.7*wn(1) & wn > 2*p.detumbleTarget;
if nnz(window) > 5
    c = polyfit(t(window), log(wn(window)), 1);
    tauSim = -1/c(1);
else
    tauSim = NaN;
end

mn = vecnorm(m, 2, 2);
satFraction = mean(any(abs(m) >= p.mtqMax - 1e-9, 2));

out.t = t;
out.rate = wn;
out.tTarget = tTarget;
out.tauSim = tauSim;
out.bdot = bd;

fprintf('tip-off rate %.2f deg/s, target %.2f deg/s, orbit period %.0f s\n', ...
    wn(1)*180/pi, p.detumbleTarget*180/pi, period);
fprintf('%-34s %12s %12s\n', '', 'hand sizing', 'simulation');
fprintf('%-34s %12.0f %12.0f\n', 'rate e-fold time [s]', bd.tauE, tauSim);
fprintf('%-34s %12.0f %12.0f\n', 'time to target rate [s]', bd.tDetumble, tTarget);
fprintf('%-34s %12.1f %12.1f\n', 'time to target rate [orbits]', bd.tDetumble/period, tTarget/period);
fprintf('B-dot gain %.3e, peak commanded dipole %.1f A m^2 (limit %.0f), saturated %.1f %% of samples\n', ...
    bd.k, max(max(abs(m))), p.mtqMax, 100*satFraction);
fprintf('final rate %.3f deg/s after %.1f orbits\n', wn(end)*180/pi, t(end)/period);

figure('Name', 'B-dot detumble');
semilogy(t/period, wn*180/pi, 'LineWidth', 1.2);
hold on;
yline(p.detumbleTarget*180/pi, '--');
grid on;
xlabel('orbits');
ylabel('|\omega| [deg/s]');
title('B-dot detumble');
if nargout == 0
    clear out
end
end
