function build_baseline_model(logStep)
% Builds models/baseline_plant.slx: J2 orbit, rigid-body attitude dynamics and the four
% disturbance torques, open loop. Integrators and signal flow are Simulink blocks; the
% force and torque models are the tested .m functions, called from MATLAB Function blocks.
if nargin < 1
    logStep = 10;
end
root = fileparts(fileparts(mfilename('fullpath')));
p = smallsat_params;
period = 2*pi/sqrt(p.mu/(p.Re + p.altitude)^3);

mdl = 'baseline_plant';
if bdIsLoaded(mdl), close_system(mdl, 0); end
oldFile = fullfile(root, 'models', [mdl '.slx']);
if exist(oldFile, 'file'), delete(oldFile); end
rehash;
new_system(mdl);
closer = onCleanup(@() close_system(mdl, 0));
set_param(mdl, 'InitFcn', 'P = smallsat_params; ic = initial_state(P);');

fcn(mdl, 'orbit_acc', {'function a = orbit_acc(r)', '%#codegen', 'a = j2_accel(r, P);', 'end'}, [250 40 350 90]);
fcn(mdl, 'disturbances', {'function [tau, parts] = disturbances(t, r, v, q)', '%#codegen', ...
    '[tau, parts] = disturbance_block(t, r, v, q, P);', 'end'}, [250 160 400 260]);
fcn(mdl, 'ang_acc', {'function wd = ang_acc(tau, w)', '%#codegen', ...
    'wd = P.inertia \ (tau - cross(w, P.inertia*w));', 'end'}, [500 160 600 220]);
fcn(mdl, 'quat_kin', {'function dq = quat_kin(q, w)', '%#codegen', 'dq = quat_rate(q, w);', 'end'}, ...
    [500 300 600 360]);

add_block('simulink/Sources/Clock', [mdl '/clock'], 'Position', [60 180 90 200]);
integ(mdl, 'vInt', 'ic.v', [420 40 450 70]);
integ(mdl, 'rInt', 'ic.r', [520 40 550 70]);
integ(mdl, 'wInt', 'ic.w', [650 160 680 190]);
integ(mdl, 'qInt', 'ic.q', [650 300 680 330]);

wire(mdl, 'rInt/1', 'orbit_acc/1');
wire(mdl, 'orbit_acc/1', 'vInt/1');
wire(mdl, 'vInt/1', 'rInt/1');
wire(mdl, 'clock/1', 'disturbances/1');
wire(mdl, 'rInt/1', 'disturbances/2');
wire(mdl, 'vInt/1', 'disturbances/3');
wire(mdl, 'qInt/1', 'disturbances/4');
wire(mdl, 'disturbances/1', 'ang_acc/1');
wire(mdl, 'wInt/1', 'ang_acc/2');
wire(mdl, 'ang_acc/1', 'wInt/1');
wire(mdl, 'qInt/1', 'quat_kin/1');
wire(mdl, 'wInt/1', 'quat_kin/2');
wire(mdl, 'quat_kin/1', 'qInt/1');

logs = {'t_log', 'clock/1'; 'r_log', 'rInt/1'; 'v_log', 'vInt/1'; 'q_log', 'qInt/1'; ...
        'w_log', 'wInt/1'; 'tau_log', 'disturbances/1'; 'parts_log', 'disturbances/2'};
for k = 1:size(logs, 1)
    name = ['to_' logs{k, 1}];
    add_block('simulink/Sinks/To Workspace', [mdl '/' name], 'VariableName', logs{k, 1}, ...
        'SaveFormat', 'Array', 'SampleTime', num2str(logStep), ...
        'Position', [800 20 + 40*k 880 40 + 40*k]);
    wire(mdl, logs{k, 2}, [name '/1']);
end

set_param(mdl, 'Solver', 'ode113', 'RelTol', '1e-10', 'AbsTol', '1e-10', ...
    'StopTime', num2str(period, 10));
modelDir = fullfile(root, 'models');
if ~exist(modelDir, 'dir'), mkdir(modelDir); end
save_system(mdl, fullfile(modelDir, [mdl '.slx']));
fprintf('saved %s\n', fullfile(modelDir, [mdl '.slx']));
end

function fcn(mdl, name, lines, pos)
blk = [mdl '/' name];
add_block('simulink/User-Defined Functions/MATLAB Function', blk, 'Position', pos);
cfg = get_param(blk, 'MATLABFunctionConfiguration');
cfg.FunctionScript = strjoin(lines, newline);
chart = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
d = Stateflow.Data(chart);
d.Name = 'P';
d.Scope = 'Parameter';
end

function integ(mdl, name, icExpr, pos)
add_block('simulink/Continuous/Integrator', [mdl '/' name], 'InitialCondition', icExpr, 'Position', pos);
end

function wire(mdl, from, to)
add_line(mdl, from, to, 'autorouting', 'on');
end
