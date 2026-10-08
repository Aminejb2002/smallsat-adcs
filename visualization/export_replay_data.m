function export_replay_data()
% Runs models/smallsat_adcs.slx twice and writes the logged signals the replay videos are drawn from
% to results/replay/ (estimate in the loop, tumble start, default noise seeds):
%   mission.csv  full three-orbit run, every 0.5 s: t, q, qref, w, h, m (dipole), tauW
%   outage.csv   the same run with a 1800 s star tracker outage from 14000 s, every 10 s: t, q, qref, qHat, sigma
% Rendering is done by visualization/render_replay.py. Run from anywhere with the project on the path.
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root, 'results', 'replay');
if ~exist(outDir, 'dir'), mkdir(outDir); end
p0 = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'tumble');
cleanup = onCleanup(@() evalin('base', 'clear adcsCase paramOverride'));

% mission, logged every 0.5 s (the model's own logging is every 10 s): only the To Workspace sample times of these
% signals are changed for this run, through the SimulationInput, the saved model is not touched
evalin('base', 'clear paramOverride');
in = Simulink.SimulationInput(mdl).setBlockParameter([mdl '/use_estimate'], 'Value', '1');
logVars = {'t_log', 'q_log', 'qref_log', 'w_log', 'h_log', 'm_log', 'tauW_log'};
load_system(mdl);
tw = find_system(mdl, 'LookUnderMasks', 'all', 'BlockType', 'ToWorkspace');
nSet = 0;
for k = 1:numel(tw)
    if any(strcmp(regexprep(get_param(tw{k}, 'VariableName'), '^out\.', ''), logVars))
        in = in.setBlockParameter(tw{k}, 'SampleTime', '0.5');
        nSet = nSet + 1;
    end
end
fprintf('fine logging: %d To Workspace blocks set to 0.5 s (expected %d)\n', nSet, numel(logVars));
out = sim(in);
t = out.t_log(:);
M = [t, log_rows(out.q_log), log_rows(out.qref_log), log_rows(out.w_log), ...
     log_rows(out.h_log), log_rows(out.m_log), log_rows(out.tauW_log)];
names = {'t', 'q0', 'q1', 'q2', 'q3', 'r0', 'r1', 'r2', 'r3', 'wx', 'wy', 'wz', ...
         'hx', 'hy', 'hz', 'mx', 'my', 'mz', 'tx', 'ty', 'tz'};
write_table(fullfile(outDir, 'mission.csv'), M, names);
fprintf('mission.csv: %d samples, %.0f s\n', size(M, 1), t(end));

% outage, model logging (10 s)
in = Simulink.SimulationInput(mdl).setBlockParameter([mdl '/use_estimate'], 'Value', '1');
assignin('base', 'paramOverride', struct('outageStart', 14000, 'outageDuration', 1800));
out = sim(in);
t = out.t_log(:);
idx = 1:round(10/p0.gyroStep):size(log_rows(out.estQ_log), 1);
qHat = log_rows(out.estQ_log);
sig = log_rows(out.estSig_log);
O = [t, log_rows(out.q_log), log_rows(out.qref_log), qHat(idx, :), sig(idx, :)];
names = {'t', 'q0', 'q1', 'q2', 'q3', 'r0', 'r1', 'r2', 'r3', 'e0', 'e1', 'e2', 'e3', 'sx', 'sy', 'sz'};
write_table(fullfile(outDir, 'outage.csv'), O, names);
fprintf('outage.csv: %d samples, outage 14000 to 15800 s\n', size(O, 1));

if ~wasLoaded, close_system(mdl, 0); end
end

function write_table(file, A, names)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', strjoin(names, ','));
fclose(fid);
dlmwrite(file, A, '-append', 'precision', '%.10g');
end
