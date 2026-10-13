function P = interference_defaults(P)
% Fills the parameters of the magnetorquer-interference study (smallsat_adcs_interference.slx) that are
% missing in P. With these defaults the interference model behaves like smallsat_adcs.slx: torquers always
% on, no leakage into the magnetometer, B-dot gain unchanged. Change them with a paramOverride struct.
%   mtqCycle    B-dot cycle length, s
%   mtqOff      torquers off at the start of each cycle for this long, s
%   mtqLeak     coupling from torquer dipole to magnetometer reading, T per A m^2 (same axis)
%   mtqLeakTau  settling time constant of the coil field after switch-off, s (0 = settled at once)
%   bdotScale   factor on the B-dot gain
d = struct('mtqCycle', 1, 'mtqOff', 0, 'mtqLeak', 0, 'mtqLeakTau', 0, 'bdotScale', 1);
for f = fieldnames(d)'
    if ~isfield(P, f{1}), P.(f{1}) = d.(f{1}); end
end
end
