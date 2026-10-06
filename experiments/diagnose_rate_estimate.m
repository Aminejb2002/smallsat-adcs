function diagnose_rate_estimate
% Three-orbit tumble run with the estimate in the loop. The estimator is replayed offline on the
% logged gyro and star tracker samples (same code as the block, see estimator_crosscheck) and the
% rate and bias it hands to the controller are compared with the true values over the last 2000 s.
% A constant rate error dw settles the pointing error at about 2*zeta/wn*dw, which is printed
% next to the measured pointing error.
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
cleanup = onCleanup(@() evalin('base', 'clear adcsCase'));
in = Simulink.SimulationInput(mdl);
in = in.setBlockParameter([mdl '/use_estimate'], 'Value', '1');
out = sim(in);

t = out.t_log(:);
w = log_rows(out.w_log);
q = log_rows(out.q_log);
qref = log_rows(out.qref_log);
wm = log_rows(out.estIn_w_log);
qm = log_rows(out.estIn_q_log);
new = log_rows(out.estIn_new_log);
tau = log_rows(out.estIn_tau_log);
hw = log_rows(out.estIn_h_log);
bTrue = log_rows(out.gyroBias_log);
n = size(wm, 1);

x = zeros(157, 1);
wHat = zeros(n, 3);
bHat = zeros(n, 3);
for k = 1:n
    x = mekf_step(x, wm(k, :)', qm(k, :)', new(k), tau(k, :)', hw(k, :)', p);
    [~, wHat(k, :), ~, bHat(k, :)] = mekf_outputs(x, wm(k, :)', qm(k, :)');
end

k = min(round(t/p.gyroStep) + 1, n);
tail = t >= t(end) - 2000;
dw = wHat(k(tail), :) - w(tail, :);
kb = min(round(t/p.gyroStep) + 1, size(bTrue, 1));
if size(bTrue, 1) == numel(t)
    kb = (1:numel(t))';
end
db = bHat(k(tail), :) - bTrue(kb(tail), :);
raw = wm(k(tail), :) - bHat(k(tail), :) - w(tail, :);

errDeg = zeros(nnz(tail), 1);
it = find(tail);
for j = 1:numel(it)
    qe = quat_mul([qref(it(j), 1); -qref(it(j), 2:4)'], q(it(j), :)');
    errDeg(j) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
end

d = 180/pi;
fprintf('last 2000 s, estimate in the loop (all rates in deg/s)\n');
fprintf('  filtered rate error  mean %s  std %s\n', mat2str(mean(dw)*d, 3), mat2str(std(dw)*d, 3));
fprintf('  unfiltered rate error mean %s  std %s\n', mat2str(mean(raw)*d, 3), mat2str(std(raw)*d, 3));
fprintf('  bias estimate error  mean %s  (true bias %s)\n', mat2str(mean(db)*d, 3), mat2str(mean(bTrue(kb(tail), :))*d, 3));
fprintf('  pointing error: measured mean %.4f deg, predicted from mean filtered rate error %.4f deg\n', ...
    mean(errDeg), 2*p.pointZeta/p.pointWn*norm(mean(dw))*d);
if ~wasLoaded
    close_system(mdl, 0);
end
end
