function diagnose_dumping
% Where the Simulink momentum dumping departs from the plain-MATLAB reference. Runs the 3000 s
% pointing case with dumping and prints, for the first logged samples and the worst one, the
% dipole of the model, of the reference, and the dipole the dump law gives from the model's own
% logged field and wheel momentum (same sample). If the model matches the third column, the
% block is right and the reference differs in timing or inputs.
p = smallsat_params;
mdl = 'smallsat_adcs';
wasLoaded = bdIsLoaded(mdl);
assignin('base', 'adcsCase', 'point');
cleanup = onCleanup(@() evalin('base', 'clear adcsCase'));
ic = initial_state_case(p, 'point');
duration = 3000;
in = Simulink.SimulationInput(mdl).setModelParameter('StopTime', num2str(duration));
in = in.setBlockParameter([mdl '/Sensors/noise_gain_B'], 'Gain', '0');
if getSimulinkBlockHandle([mdl '/use_estimate']) > 0
    in = in.setBlockParameter([mdl '/use_estimate'], 'Value', '0');
end
out = sim(in);
ref = reference_run(duration, 10, 1e-12, ic, struct('pointing', true));
t = out.t_log(:);
B = log_rows(out.B_log);
h = log_rows(out.h_log);
m = log_rows(out.m_log);
mFrom = zeros(size(m));
for k = 1:numel(t)
    mFrom(k, :) = dump_law(B(k, :)', h(k, :)', p)';
end
dRef = vecnorm(m - ref.m, 2, 2);
dSelf = vecnorm(m - mFrom, 2, 2);
dRefSelf = vecnorm(ref.m - mFrom, 2, 2);
fprintf('max |model m - reference m|              %.3e A m^2\n', max(dRef));
fprintf('max |model m - law(model B, model h)|    %.3e A m^2\n', max(dSelf));
fprintf('max |reference m - law(model B, model h)| %.3e A m^2\n', max(dRefSelf));
fprintf('\n%6s %26s %26s %26s\n', 't [s]', 'model m', 'reference m', 'law(model B, h)');
for k = [1:6 find(dRef == max(dRef), 1)]
    fprintf('%6.0f  %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f   %7.3f %7.3f %7.3f\n', t(k), m(k, :), ref.m(k, :), mFrom(k, :));
end
fprintf('\n|h| model %.4e  reference %.4e at t = %.0f s\n', norm(h(end, :)), norm(ref.h(end, :)), t(end));
if ~wasLoaded
    close_system(mdl, 0);
end
end
