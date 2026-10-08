function build_adcs_model(logStep, force)
% Builds the scaffold models/smallsat_adcs.slx and a backup copy smallsat_adcs_scaffold.slx.
% Plant (orbit, attitude, three reaction wheels, disturbances), magnetometer, guidance and
% actuators are built here. The Controller subsystem holds the B-dot loop and placeholders for
% the pointing law and the mode switch, which are built by hand in Simulink.
% Refuses to overwrite an existing model unless force is true. This only documents how the first
% scaffold was made: the model in the repository was extended by hand afterwards (controller, estimator,
% sensors, logging, masks) and no longer matches what this script generates.
if nargin < 1
    logStep = 10;
end
if nargin < 2
    force = false;
end
root = fileparts(fileparts(mfilename('fullpath')));
p = smallsat_params;
period = 2*pi*sqrt((p.Re + p.altitude)^3/p.mu);

mdl = 'smallsat_adcs';
modelDir = fullfile(root, 'models');
modelFile = fullfile(modelDir, [mdl '.slx']);
if exist(modelFile, 'file') && ~force
    error('build_adcs_model:exists', ...
        '%s exists and may hold hand-built work. Call build_adcs_model(10, true) to overwrite it.', modelFile);
end
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(modelFile, 'file'), delete(modelFile); end
rehash;
new_system(mdl);
closer = onCleanup(@() close_system(mdl, 0));
set_param(mdl, 'InitFcn', ['P = smallsat_params; ' ...
    'if ~exist(''adcsCase'', ''var''), adcsCase = ''tumble''; end; ' ...
    'ic = initial_state_case(P, adcsCase); bdot = size_bdot(P); pointing = size_pointing(P);']);

build_plant([mdl '/Plant'], [380 60 560 340]);
build_sensors([mdl '/Sensors'], [680 60 800 110]);
build_guidance([mdl '/Guidance'], [380 400 500 450]);
build_controller([mdl '/Controller'], [680 180 820 300]);
build_actuators([mdl '/Actuators'], [680 340 820 420]);

wire(mdl, 'Plant/5', 'Sensors/1');
wire(mdl, 'Plant/1', 'Guidance/1');
wire(mdl, 'Plant/2', 'Guidance/2');
wire(mdl, 'Sensors/1', 'Controller/1');
wire(mdl, 'Plant/3', 'Controller/2');
wire(mdl, 'Plant/4', 'Controller/3');
wire(mdl, 'Guidance/1', 'Controller/4');
wire(mdl, 'Guidance/2', 'Controller/5');
wire(mdl, 'Controller/1', 'Actuators/1');
wire(mdl, 'Controller/2', 'Actuators/2');
wire(mdl, 'Plant/8', 'Actuators/3');
wire(mdl, 'Actuators/1', 'Plant/1');
wire(mdl, 'Actuators/2', 'Plant/2');

add_block('simulink/Sources/Clock', [mdl '/clock'], 'Position', [60 40 90 60]);
logs = {'t_log', 'clock/1'; 'r_log', 'Plant/1'; 'v_log', 'Plant/2'; 'q_log', 'Plant/3'; ...
        'w_log', 'Plant/4'; 'B_log', 'Plant/5'; 'parts_log', 'Plant/6'; 'tau_mt_log', 'Plant/7'; ...
        'h_log', 'Plant/8'; 'Bm_log', 'Sensors/1'; 'qref_log', 'Guidance/1'; ...
        'm_log', 'Actuators/1'; 'tauW_log', 'Actuators/2'};
for k = 1:size(logs, 1)
    name = ['to_' logs{k, 1}];
    add_block('simulink/Sinks/To Workspace', [mdl '/' name], 'VariableName', logs{k, 1}, ...
        'SaveFormat', 'Array', 'SampleTime', num2str(logStep), ...
        'Position', [920 20 + 36*k 1000 40 + 36*k]);
    wire(mdl, logs{k, 2}, [name '/1']);
end

set_param(mdl, 'Solver', 'ode113', 'RelTol', '1e-10', 'AbsTol', '1e-10', ...
    'StopTime', num2str(3*round(period)));
if ~exist(modelDir, 'dir'), mkdir(modelDir); end
save_system(mdl, modelFile);
copyfile(modelFile, fullfile(modelDir, [mdl '_scaffold.slx']));
fprintf('saved %s and a backup copy %s_scaffold.slx\n', modelFile, mdl);
end

function build_plant(sys, pos)
newSystem(sys, pos);
add_block('simulink/Sources/In1', [sys '/m'], 'Position', [30 140 60 154]);
add_block('simulink/Sources/In1', [sys '/tauW'], 'Position', [30 200 60 214]);
add_block('simulink/Sources/Clock', [sys '/clock'], 'Position', [30 40 60 60]);

fcn(sys, 'orbit_acc', {'function a = orbit_acc(r)', '%#codegen', 'a = j2_accel(r, P);', 'end'}, ...
    [250 20 350 70], true);
fcn(sys, 'disturbances', {'function [tau, parts] = disturbances(t, r, v, q)', '%#codegen', ...
    '[tau, parts] = disturbance_block(t, r, v, q, P);', 'end'}, [250 110 400 210], true);
fcn(sys, 'mag_body', {'function B = mag_body(t, r, q)', '%#codegen', ...
    'B = quat_to_dcm(q)*magnetic_field_eci(r, t, P);', 'end'}, [250 250 350 310], true);
fcn(sys, 'mt_torque', {'function tau = mt_torque(m, B)', '%#codegen', 'tau = cross(m, B);', 'end'}, ...
    [450 250 540 300], false);
fcn(sys, 'ang_acc', {'function wd = ang_acc(tau, w, h, tauW)', '%#codegen', ...
    'wd = angular_accel(w, h, tau, tauW, P.inertia);', 'end'}, [640 130 740 210], true);
fcn(sys, 'quat_kin', {'function dq = quat_kin(q, w)', '%#codegen', 'dq = quat_rate(q, w);', 'end'}, ...
    [640 350 740 410], false);
add_block('simulink/Math Operations/Sum', [sys '/tau_sum'], 'Inputs', '++', 'Position', [500 140 520 180]);

integ(sys, 'vInt', 'ic.v', [420 20 450 50]);
integ(sys, 'rInt', 'ic.r', [520 20 550 50]);
integ(sys, 'wInt', 'ic.w', [800 130 830 160]);
integ(sys, 'qInt', 'ic.q', [800 350 830 380]);
integ(sys, 'hInt', 'ic.h', [800 230 830 260]);

outs = {'r', 'v', 'q', 'w', 'B', 'parts', 'tau_mt', 'h'};
for k = 1:numel(outs)
    add_block('simulink/Sinks/Out1', [sys '/' outs{k}], 'Position', [980 20 + 50*k 1010 34 + 50*k]);
end

wire(sys, 'rInt/1', 'orbit_acc/1');
wire(sys, 'orbit_acc/1', 'vInt/1');
wire(sys, 'vInt/1', 'rInt/1');
wire(sys, 'clock/1', 'disturbances/1');
wire(sys, 'rInt/1', 'disturbances/2');
wire(sys, 'vInt/1', 'disturbances/3');
wire(sys, 'qInt/1', 'disturbances/4');
wire(sys, 'clock/1', 'mag_body/1');
wire(sys, 'rInt/1', 'mag_body/2');
wire(sys, 'qInt/1', 'mag_body/3');
wire(sys, 'm/1', 'mt_torque/1');
wire(sys, 'mag_body/1', 'mt_torque/2');
wire(sys, 'disturbances/1', 'tau_sum/1');
wire(sys, 'mt_torque/1', 'tau_sum/2');
wire(sys, 'tau_sum/1', 'ang_acc/1');
wire(sys, 'wInt/1', 'ang_acc/2');
wire(sys, 'hInt/1', 'ang_acc/3');
wire(sys, 'tauW/1', 'ang_acc/4');
wire(sys, 'ang_acc/1', 'wInt/1');
wire(sys, 'tauW/1', 'hInt/1');
wire(sys, 'qInt/1', 'quat_kin/1');
wire(sys, 'wInt/1', 'quat_kin/2');
wire(sys, 'quat_kin/1', 'qInt/1');
wire(sys, 'rInt/1', 'r/1');
wire(sys, 'vInt/1', 'v/1');
wire(sys, 'qInt/1', 'q/1');
wire(sys, 'wInt/1', 'w/1');
wire(sys, 'mag_body/1', 'B/1');
wire(sys, 'disturbances/2', 'parts/1');
wire(sys, 'mt_torque/1', 'tau_mt/1');
wire(sys, 'hInt/1', 'h/1');
end

function build_sensors(sys, pos)
newSystem(sys, pos);
add_block('simulink/Sources/In1', [sys '/B'], 'Position', [30 40 60 54]);
add_block('simulink/Discrete/Zero-Order Hold', [sys '/hold'], 'SampleTime', 'P.sensorStep', ...
    'Position', [120 30 160 60]);
add_block('simulink/Sources/Random Number', [sys '/noise'], 'Mean', '0', 'Variance', '1', ...
    'Seed', '[11 22 33]', 'SampleTime', 'P.sensorStep', 'Position', [60 120 110 150]);
add_block('simulink/Math Operations/Reshape', [sys '/to_column'], 'OutputDimensionality', ...
    'Customize', 'OutputDimensions', '[3 1]', 'Position', [140 120 180 150]);
add_block('simulink/Math Operations/Gain', [sys '/noise_gain'], 'Gain', 'P.magSigma', ...
    'Position', [210 120 250 150]);
add_block('simulink/Math Operations/Sum', [sys '/add_noise'], 'Inputs', '++', 'Position', [300 40 320 80]);
add_block('simulink/Sinks/Out1', [sys '/Bmeas'], 'Position', [400 52 430 66]);
wire(sys, 'B/1', 'hold/1');
wire(sys, 'hold/1', 'add_noise/1');
wire(sys, 'noise/1', 'to_column/1');
wire(sys, 'to_column/1', 'noise_gain/1');
wire(sys, 'noise_gain/1', 'add_noise/2');
wire(sys, 'add_noise/1', 'Bmeas/1');
end

function build_guidance(sys, pos)
newSystem(sys, pos);
add_block('simulink/Sources/In1', [sys '/r'], 'Position', [30 40 60 54]);
add_block('simulink/Sources/In1', [sys '/v'], 'Position', [30 100 60 114]);
fcn(sys, 'nadir', {'function [qref, wref] = nadir(r, v)', '%#codegen', ...
    '[qref, wref] = nadir_reference(r, v);', 'end'}, [140 40 240 120], false);
add_block('simulink/Sinks/Out1', [sys '/qref'], 'Position', [340 52 370 66]);
add_block('simulink/Sinks/Out1', [sys '/wref'], 'Position', [340 92 370 106]);
wire(sys, 'r/1', 'nadir/1');
wire(sys, 'v/1', 'nadir/2');
wire(sys, 'nadir/1', 'qref/1');
wire(sys, 'nadir/2', 'wref/1');
end

function build_controller(sys, pos)
newSystem(sys, pos);
inputs = {'Bmeas', 'q', 'w', 'qref', 'wref'};
for k = 1:numel(inputs)
    add_block('simulink/Sources/In1', [sys '/' inputs{k}], 'Position', [30 40 + 60*k 60 54 + 60*k]);
end
for k = 2:numel(inputs)
    add_block('simulink/Sinks/Terminator', [sys '/unused_' inputs{k}], ...
        'Position', [110 40 + 60*k 130 60 + 60*k]);
    wire(sys, [inputs{k} '/1'], ['unused_' inputs{k} '/1']);
end
add_block('simulink/Discrete/Unit Delay', [sys '/previous'], ...
    'InitialCondition', 'quat_to_dcm(ic.q)*magnetic_field_eci(ic.r, 0, P)', ...
    'SampleTime', 'P.sensorStep', 'Position', [220 60 260 90]);
add_block('simulink/Math Operations/Sum', [sys '/difference'], 'Inputs', '+-', 'Position', [320 20 340 60]);
add_block('simulink/Math Operations/Gain', [sys '/per_step'], 'Gain', '1/P.sensorStep', ...
    'Position', [390 25 440 55]);
add_block('simulink/Math Operations/Gain', [sys '/bdot_gain'], 'Gain', '-bdot.k', ...
    'Position', [490 25 540 55]);
add_block('simulink/Sinks/Out1', [sys '/m_cmd'], 'Position', [620 32 650 46]);
add_block('simulink/Sources/Constant', [sys '/tauW_placeholder'], 'Value', 'zeros(3, 1)', ...
    'VectorParams1D', 'off', 'Position', [490 160 540 190]);
add_block('simulink/Sinks/Out1', [sys '/tauW_cmd'], 'Position', [620 168 650 182]);
wire(sys, 'Bmeas/1', 'difference/1');
wire(sys, 'Bmeas/1', 'previous/1');
wire(sys, 'previous/1', 'difference/2');
wire(sys, 'difference/1', 'per_step/1');
wire(sys, 'per_step/1', 'bdot_gain/1');
wire(sys, 'bdot_gain/1', 'm_cmd/1');
wire(sys, 'tauW_placeholder/1', 'tauW_cmd/1');
end

function build_actuators(sys, pos)
newSystem(sys, pos);
add_block('simulink/Sources/In1', [sys '/m_cmd'], 'Position', [30 40 60 54]);
add_block('simulink/Sources/In1', [sys '/tauW_cmd'], 'Position', [30 120 60 134]);
add_block('simulink/Sources/In1', [sys '/h'], 'Position', [30 160 60 174]);
add_block('simulink/Discontinuities/Saturation', [sys '/torquer_limit'], 'UpperLimit', 'P.mtqMax', ...
    'LowerLimit', '-P.mtqMax', 'Position', [120 30 170 65]);
fcn(sys, 'wheel_motor_block', {'function tauW = wheel_motor_block(cmd, h)', '%#codegen', ...
    'tauW = wheel_motor(cmd, h, P);', 'end'}, [120 110 240 180], true);
add_block('simulink/Sinks/Out1', [sys '/m'], 'Position', [320 40 350 54]);
add_block('simulink/Sinks/Out1', [sys '/tauW'], 'Position', [320 138 350 152]);
wire(sys, 'm_cmd/1', 'torquer_limit/1');
wire(sys, 'torquer_limit/1', 'm/1');
wire(sys, 'tauW_cmd/1', 'wheel_motor_block/1');
wire(sys, 'h/1', 'wheel_motor_block/2');
wire(sys, 'wheel_motor_block/1', 'tauW/1');
end

function newSystem(sys, pos)
add_block('built-in/Subsystem', sys, 'Position', pos);
Simulink.SubSystem.deleteContents(sys);
end

function fcn(sys, name, lines, pos, needsP)
blk = [sys '/' name];
add_block('simulink/User-Defined Functions/MATLAB Function', blk, 'Position', pos);
cfg = get_param(blk, 'MATLABFunctionConfiguration');
cfg.FunctionScript = strjoin(lines, newline);
if needsP
    chart = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
    d = Stateflow.Data(chart);
    d.Name = 'P';
    d.Scope = 'Parameter';
end
end

function integ(sys, name, icExpr, pos)
add_block('simulink/Continuous/Integrator', [sys '/' name], 'InitialCondition', icExpr, 'Position', pos);
end

function wire(sys, from, to)
add_line(sys, from, to, 'autorouting', 'on');
end
