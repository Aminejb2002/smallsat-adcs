function probe_matlab_function
% Checks that a build script can create MATLAB Function blocks, give them the
% parameter struct, and close an integrator loop; writes results/probe_matlab_function.txt.

root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
reportFile = fullfile(outDir, 'probe_matlab_function.txt');
fid = fopen(reportFile, 'w');
closer = onCleanup(@() fclose(fid));

p = smallsat_params;
assignin('base', 'P', p);
a0 = p.Re + p.altitude;
[r0, v0] = circular_state(a0, sso_inclination(a0, p), 0.3, 0, p);

mdl = 'probe_fcn_model';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
c1 = onCleanup(@() close_system(mdl, 0));
build_fcn_block(mdl, fid);
add_block('simulink/Sources/Constant', [mdl '/r'], 'Value', mat2str(r0));
add_block('simulink/Sinks/To Workspace', [mdl '/out'], 'VariableName', 'a_probe', 'SaveFormat', 'Array');
add_line(mdl, 'r/1', 'fcn/1');
add_line(mdl, 'fcn/1', 'out/1');
set_param(mdl, 'Solver', 'FixedStepDiscrete', 'FixedStep', '1', 'StopTime', '0');
try
    simOut = sim(mdl);
    rows = as_rows(simOut.a_probe);
    got = rows(1, :)';
    expected = j2_accel(r0, p);
    relErr = max(abs(got - expected))/max(abs(expected));
    emit(fid, '%-8s block output vs j2_accel, relative difference %.2e', status(relErr < 1e-12), relErr);
catch e
    emit(fid, 'FAILED   simulation of single block: %s', e.message);
end

mdl2 = 'probe_orbit_model';
if bdIsLoaded(mdl2), close_system(mdl2, 0); end
new_system(mdl2);
c2 = onCleanup(@() close_system(mdl2, 0));
build_fcn_block(mdl2, fid);
add_block('simulink/Continuous/Integrator', [mdl2 '/vInt'], 'InitialCondition', mat2str(v0));
add_block('simulink/Continuous/Integrator', [mdl2 '/rInt'], 'InitialCondition', mat2str(r0));
add_block('simulink/Sinks/To Workspace', [mdl2 '/out'], 'VariableName', 'r_probe', 'SaveFormat', 'Array');
add_line(mdl2, 'fcn/1', 'vInt/1');
add_line(mdl2, 'vInt/1', 'rInt/1');
add_line(mdl2, 'rInt/1', 'fcn/1');
add_line(mdl2, 'rInt/1', 'out/1');
set_param(mdl2, 'Solver', 'ode113', 'RelTol', '1e-10', 'AbsTol', '1e-6', 'StopTime', '600');
try
    simOut = sim(mdl2);
    rows = as_rows(simOut.r_probe);
    rSim = rows(end, :)';
    opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-6);
    [~, x] = ode113(@(t, x) [x(4:6); j2_accel(x(1:3), p)], [0 600], [r0; v0], opts);
    dr = norm(rSim - x(end, 1:3)');
    emit(fid, '%-8s integrator loop vs ode113 after 600 s, position difference %.3e m', status(dr < 0.01), dr);
catch e
    emit(fid, 'FAILED   simulation of integrator loop: %s', e.message);
end

emit(fid, '\nReport written to %s', reportFile);
end

function build_fcn_block(mdl, fid)
blk = [mdl '/fcn'];
add_block('simulink/User-Defined Functions/MATLAB Function', blk);
code = sprintf('function a = fcn(r)\n%%#codegen\na = j2_accel(r, P);\nend\n');
scriptSet = false;
try
    cfg = get_param(blk, 'MATLABFunctionConfiguration');
    cfg.FunctionScript = code;
    scriptSet = true;
    emit(fid, 'ok       script set through MATLABFunctionConfiguration');
catch e
    emit(fid, 'FAILED   MATLABFunctionConfiguration: %s', e.message);
end
chart = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', blk);
if ~scriptSet
    try
        chart.Script = code;
        emit(fid, 'ok       script set through Stateflow EMChart');
    catch e
        emit(fid, 'FAILED   EMChart script: %s', e.message);
    end
end
try
    d = Stateflow.Data(chart);
    d.Name = 'P';
    d.Scope = 'Parameter';
    emit(fid, 'ok       parameter data P added');
catch e
    emit(fid, 'FAILED   parameter data: %s', e.message);
end
end

function rows = as_rows(a)
if ndims(a) == 3 || (size(a, 1) == 3 && size(a, 2) == 1)
    rows = reshape(a, 3, [])';
else
    rows = a;
end
end

function s = status(ok)
if ok
    s = 'ok';
else
    s = 'MISMATCH';
end
end

function emit(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});
fprintf('%s\n', msg);
fprintf(fid, '%s\n', msg);
end
