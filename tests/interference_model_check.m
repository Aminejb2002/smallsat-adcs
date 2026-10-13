function ok = interference_model_check(stopTime)
% With default parameters (torquers always on, no leakage, gain unchanged) the interference copy of the model
% must give the same results as smallsat_adcs.slx. Runs both from the tumble case and compares the dipole,
% body rate, magnetometer sample and attitude logs. Build the copy first with build_interference_model.
if nargin < 1, stopTime = 1200; end
names = {'smallsat_adcs', 'smallsat_adcs_interference'};
logs = {'m_log', 'w_log', 'Bm_log', 'q_log'};
cleanup = onCleanup(@() evalin('base', 'clear adcsCase paramOverride'));
assignin('base', 'adcsCase', 'tumble');
evalin('base', 'if exist(''paramOverride'', ''var''), clear paramOverride; end');
res = cell(1, 2);
for i = 1:2
    wasLoaded = bdIsLoaded(names{i});
    load_system(names{i});
    out = sim(names{i}, 'StopTime', num2str(stopTime));
    for j = 1:numel(logs)
        res{i}.(logs{j}) = log_rows(out.(logs{j}));
    end
    if ~wasLoaded, close_system(names{i}, 0); end
end
ok = true;
for j = 1:numel(logs)
    a = res{1}.(logs{j}); b = res{2}.(logs{j});
    if ~isequal(size(a), size(b))
        fprintf('%-7s size differs: %s vs %s\n', logs{j}, mat2str(size(a)), mat2str(size(b)));
        ok = false; continue
    end
    d = max(abs(a(:) - b(:)));
    scale = max(abs(a(:)));
    fprintf('%-7s max difference %.3e (largest value %.3e)\n', logs{j}, d, scale);
    ok = ok && d <= 1e-9*max(scale, 1e-12);
end
if ok, fprintf('interference_model_check: copy matches the original with default parameters\n');
else, fprintf('interference_model_check: DIFFERENCE FOUND\n'); end
end
