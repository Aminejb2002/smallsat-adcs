function build_interference_model(force)
% Builds models/smallsat_adcs_interference.slx as a COPY of models/smallsat_adcs.slx with
% magnetorquer-magnetometer interference and a sensing / actuation cycle for B-dot. The original model
% file is only read, never saved. Changes in the copy:
%   Controller   the B-dot path (difference, per_step, bdot_gain) is replaced by the MATLAB Function block
%                bdot_cycle (src/control/bdot_cycle.m); bdot_gain is left in place but terminated.
%                New output 3, leak_dipole.
%   Sensors      new input 4, leak_dipole, scaled by P.mtqLeak and added to the magnetometer reading.
%   top level    leak_dipole goes through a one-sample delay (breaks the algebraic loop) into Sensors.
%   InitFcn      P = interference_defaults(P) after P = smallsat_params.
% With the defaults of interference_defaults the copy should give the same results as the original; check
% that with tests/interference_model_check.m before using it.
if nargin < 1, force = false; end
root = fileparts(fileparts(mfilename('fullpath')));
src = 'smallsat_adcs';
dst = 'smallsat_adcs_interference';
dstFile = fullfile(root, 'models', [dst '.slx']);
if exist(dstFile, 'file') && ~force
    error('build_interference_model:exists', '%s exists. Call build_interference_model(true) to rebuild it.', dstFile);
end
if bdIsLoaded(dst), close_system(dst, 0); end
if bdIsLoaded(src), close_system(src, 0); end
load_system(fullfile(root, 'models', [src '.slx']));
save_system(src, dstFile);          % from here on the loaded model is the copy; the original file is untouched
mdl = dst;

% InitFcn
init = get_param(mdl, 'InitFcn');
assert(contains(init, 'P = smallsat_params;'), 'InitFcn has an unexpected form');
set_param(mdl, 'InitFcn', strrep(init, 'P = smallsat_params;', 'P = smallsat_params; P = interference_defaults(P);'));

% Controller
ctl = [mdl '/Controller'];
lh = get_param([ctl '/bdot_gain'], 'LineHandles');
delete_line(lh.Outport);
add_block('simulink/Sinks/Terminator', [ctl '/bdot_gain_unused'], 'Position', [1100 560 1120 580]);
add_line(ctl, 'bdot_gain/1', 'bdot_gain_unused/1');
add_block('simulink/Sources/Constant', [ctl '/bdot_k'], 'Value', 'bdot.k', 'Position', [900 600 960 620]);
txt = fileread(fullfile(root, 'src', 'control', 'bdot_cycle.m'));
txt = strrep(txt, 'k, P)', 'k)');    % in the block, P is parameter data and not an argument
blk = [ctl '/bdot_cycle'];
add_block('simulink/User-Defined Functions/MATLAB Function', blk, 'Position', [1000 540 1080 620]);
cfg = get_param(blk, 'MATLABFunctionConfiguration');
cfg.FunctionScript = txt;
chart = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
d = Stateflow.Data(chart);
d.Name = 'P';
d.Scope = 'Parameter';
add_line(ctl, 'Bmeas/1', 'bdot_cycle/1', 'autorouting', 'on');
add_line(ctl, 'previous/1', 'bdot_cycle/2', 'autorouting', 'on');
add_line(ctl, 'bdot_k/1', 'bdot_cycle/3', 'autorouting', 'on');
add_line(ctl, 'bdot_cycle/1', 'Switch1/3', 'autorouting', 'on');
add_block('simulink/Sinks/Out1', [ctl '/leak_dipole'], 'Port', '3', 'Position', [1150 600 1180 620]);
add_line(ctl, 'bdot_cycle/2', 'leak_dipole/1', 'autorouting', 'on');

% Sensors
sen = [mdl '/Sensors'];
add_block('simulink/Sources/In1', [sen '/leak_dipole'], 'Port', '4', 'Position', [20 700 50 714]);
add_block('simulink/Math Operations/Gain', [sen '/leak_c'], 'Gain', 'P.mtqLeak', 'Position', [100 695 150 720]);
add_block('simulink/Math Operations/Sum', [sen '/add_leak'], 'Inputs', '++', 'Position', [900 40 920 80]);
lh = get_param([sen '/add_noise'], 'LineHandles');
delete_line(lh.Outport);
add_line(sen, 'add_noise/1', 'add_leak/1', 'autorouting', 'on');
add_line(sen, 'leak_dipole/1', 'leak_c/1', 'autorouting', 'on');
add_line(sen, 'leak_c/1', 'add_leak/2', 'autorouting', 'on');
add_line(sen, 'add_leak/1', 'Bmeas/1', 'autorouting', 'on');

% top level
add_block('simulink/Discrete/Unit Delay', [mdl '/leak_delay'], 'SampleTime', 'P.sensorStep', ...
    'InitialCondition', 'zeros(3,1)', 'Position', [560 600 600 630]);
add_line(mdl, 'Controller/3', 'leak_delay/1', 'autorouting', 'on');
add_line(mdl, 'leak_delay/1', 'Sensors/4', 'autorouting', 'on');

save_system(mdl);
close_system(mdl, 0);
fprintf('built %s (original %s.slx was not saved)\n', dstFile, src);
end
