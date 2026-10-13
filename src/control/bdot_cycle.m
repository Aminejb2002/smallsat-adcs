function [m, leak] = bdot_cycle(Bm, Bprev1, k, P)
% B-dot with an optional sensing / actuation cycle, called once per sensor sample (1 s).
% Used as a MATLAB Function block in smallsat_adcs_interference.slx (the block has no P argument, P is its
% parameter data) and by tests/bdot_cycle_check.m.
%   Bm      magnetometer sample (T), body frame, may contain the leakage
%   Bprev1  sample one step earlier from the model's UnitDelay; at the first call it holds the true start field,
%           which makes the first dipole identical to the original model
%   k       B-dot gain from size_bdot (A m^2 s/T)
% P fields (defaults in interference_defaults):
%   mtqCycle  cycle length, s (1 = torquers always on, as in the original model)
%   mtqOff    torquers off for the first mtqOff s of each cycle; the field is read at the end of that window
%   mtqLeakTau  settling time constant of the coil field after switch-off, s (0 = settled at once)
%   bdotScale factor on the gain
% m     dipole to apply (zero during the off window, the held command otherwise)
% leak  dipole that the magnetometer 'feels' at the next read: the last commanded value, decayed by
%       exp(-mtqOff/mtqLeakTau) over the off window (1 without an off window). Multiplied by the coupling
%       P.mtqLeak in the Sensors subsystem.
persistent Bold mHold n started
if isempty(started)
    started = true;
    Bold = Bprev1;
    mHold = zeros(3, 1);
    n = 0;
end
N = max(1, round(P.mtqCycle/P.sensorStep));
nOff = round(P.mtqOff/P.sensorStep);
phase = mod(n, N);
if phase == mod(nOff, N)
    mHold = -k*P.bdotScale*(Bm - Bold)/(N*P.sensorStep);
    mHold = min(max(mHold, -P.mtqMax), P.mtqMax);
    Bold = Bm;
end
if phase < nOff
    m = zeros(3, 1);
else
    m = mHold;
end
if nOff == 0
    f = 1;
elseif P.mtqLeakTau > 0
    f = exp(-P.mtqOff/P.mtqLeakTau);
else
    f = 0;
end
leak = f*mHold;
n = n + 1;
end
