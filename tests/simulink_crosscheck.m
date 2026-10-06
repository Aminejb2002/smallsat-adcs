function simulink_crosscheck(withPointing)
% Runs models/smallsat_adcs.slx and compares it with the plain-MATLAB reference propagated at
% tighter tolerance. Cases: open loop (torquers off) and B-dot loop for one orbit; with
% simulink_crosscheck(true) also the wheel pointing capture. Magnetometer noise is off.
% Thresholds were fixed before the first run. The model is loaded, never rebuilt.
if nargin < 1
    withPointing = false;
end
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
fid = fopen(fullfile(outDir, 'simulink_crosscheck.txt'), 'w');
closer = onCleanup(@() fclose(fid));

mdl = 'smallsat_adcs';
if ~exist(fullfile(root, 'models', [mdl '.slx']), 'file')
    error('simulink_crosscheck:nomodel', 'models/%s.slx not found, run build_adcs_model(10, true) first', mdl);
end
wasLoaded = bdIsLoaded(mdl);
p = smallsat_params;
logStep = 10;
bd = size_bdot(p);
cleanup = onCleanup(@() evalin('base', 'clear adcsCase'));

base = Simulink.SimulationInput(mdl);
base = base.setBlockParameter([mdl '/Sensors/noise_gain_B'], 'Gain', '0');
if getSimulinkBlockHandle([mdl '/use_estimate']) > 0
    base = base.setBlockParameter([mdl '/use_estimate'], 'Value', '0');
end

fails = 0;
total = 0;

assignin('base', 'adcsCase', 'tumble');
ic = initial_state_case(p, 'tumble');
duration = 5800;
in = base.setModelParameter('StopTime', num2str(duration));
openLoop = in.setBlockParameter([mdl '/Actuators/torquer_limit'], 'UpperLimit', '0');
openLoop = openLoop.setBlockParameter([mdl '/Actuators/torquer_limit'], 'LowerLimit', '0');
[f, n] = compare_case(sim(openLoop), reference_run(duration, logStep, 1e-12, ic, ...
    struct('gain', 0, 'mMax', 0)), 'open loop', fid);
fails = fails + f;
total = total + n;

[f, n] = compare_case(sim(in), reference_run(duration, logStep, 1e-12, ic, ...
    struct('gain', bd.k, 'mMax', p.mtqMax)), 'B-dot loop', fid);
fails = fails + f;
total = total + n;

if withPointing
    assignin('base', 'adcsCase', 'point');
    ic = initial_state_case(p, 'point');
    duration = 1500;
    in = base.setModelParameter('StopTime', num2str(duration));
    [f, n] = compare_case(sim(in), reference_run(duration, logStep, 1e-12, ic, ...
        struct('pointing', true)), 'wheel pointing', fid);
    fails = fails + f;
    total = total + n;
end

emit(fid, '\n%d of %d checks passed', total - fails, total);
if ~wasLoaded
    close_system(mdl, 0);
end
if fails > 0
    error('simulink_crosscheck:failed', '%d checks failed', fails);
end
end

function [fails, total] = compare_case(simOut, ref, label, fid)
t = simOut.t_log(:);
if numel(t) ~= numel(ref.t) || max(abs(t - ref.t)) > 1e-6
    error('simulink_crosscheck:grid', 'time grids differ: %d samples vs %d', numel(t), numel(ref.t));
end
r = log_rows(simOut.r_log);
v = log_rows(simOut.v_log);
q = log_rows(simOut.q_log);
w = log_rows(simOut.w_log);
h = log_rows(simOut.h_log);
B = log_rows(simOut.B_log);
m = log_rows(simOut.m_log);
tauW = log_rows(simOut.tauW_log);
parts = log_rows(simOut.parts_log);

sgn = sign(sum(q.*ref.q, 2));
checks = {
    'position',          max(vecnorm(r - ref.r, 2, 2)),                  0.05,   'm'
    'velocity',          max(vecnorm(v - ref.v, 2, 2)),                  1e-4,   'm/s'
    'attitude',          max(2*vecnorm(q.*sgn - ref.q, 2, 2)),           1e-6,   'rad'
    'body rate',         max(vecnorm(w - ref.w, 2, 2)),                  1e-9,   'rad/s'
    'wheel momentum',    max(vecnorm(h - ref.h, 2, 2)),                  1e-8,   'N m s'
    'wheel torque',      max(vecnorm(tauW - ref.tauW, 2, 2)),            1e-8,   'N m'
    'body field',        max(vecnorm(B - ref.B, 2, 2)),                  1e-11,  'T'
    'torquer dipole',    max(vecnorm(m - ref.m, 2, 2)),                  1e-6,   'A m^2'
    'torque sources',    max(abs(parts - ref.parts), [], 'all'),         1e-10,  'N m'
    'quaternion norm',   max(abs(vecnorm(q, 2, 2) - 1)),                 1e-8,   ''
    };
fails = 0;
total = size(checks, 1);
emit(fid, '\n%s: %d samples over %.0f s, Simulink ode113 1e-10 vs MATLAB ode113 1e-12', ...
    label, numel(t), t(end));
for k = 1:total
    ok = checks{k, 2} < checks{k, 3};
    if ok
        status = 'PASS';
    else
        status = 'FAIL';
        fails = fails + 1;
    end
    emit(fid, '%s  %-16s max difference %.3e %s (limit %.0e)', status, checks{k, 1}, checks{k, 2}, ...
        checks{k, 4}, checks{k, 3});
end
end

function emit(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('%s\n', msg);
fprintf(fid, '%s\n', msg);
end
