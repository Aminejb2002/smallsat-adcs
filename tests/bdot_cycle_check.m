function ok = bdot_cycle_check
% Checks the logic of bdot_cycle (the block used in smallsat_adcs_interference.slx) against the sampled
% B-dot loop in experiments/run_detumble_interference.m, on a synthetic field that does not depend on the
% dipole, so only the cycle, hold and leakage bookkeeping is compared. No Simulink needed.
P = struct('sensorStep', 1, 'mtqMax', 30, 'mtqCycle', 1, 'mtqOff', 0, 'mtqLeak', 0, 'mtqLeakTau', 0, 'bdotScale', 1);
k = 8.3e6;
cases = [1 0 0 0; 1 0 1e-6 0; 4 1 1e-6 0.5; 4 1 1e-6 0; 10 2 1e-7 0.3; 6 2 1e-6 1.0];   % cycle, off, c, tau
nSteps = 400;
t = (0:nSteps-1)';
Btrue = 3e-5*[sin(0.011*t), cos(0.017*t), 0.5*sin(0.023*t + 1)];
ok = true;
for i = 1:size(cases, 1)
    P.mtqCycle = cases(i, 1); P.mtqOff = cases(i, 2); c = cases(i, 3); P.mtqLeakTau = cases(i, 4);
    P.mtqLeak = c;
    % block, wired as in the model: leak goes through a one-sample delay into the magnetometer
    clear bdot_cycle
    mBlock = zeros(nSteps, 3); leakDelayed = zeros(3, 1); Bprev = Btrue(1, :)';
    for n = 1:nSteps
        Bm = Btrue(n, :)' + P.mtqLeak*leakDelayed;
        [m, leak] = bdot_cycle(Bm, Bprev, k, P);
        mBlock(n, :) = m';
        leakDelayed = leak;
        Bprev = Bm;
    end
    % reference: the loop of run_detumble_interference
    N = P.mtqCycle; nOff = P.mtqOff; tauS = P.mtqLeakTau;
    mRef = zeros(nSteps, 3); mHold = zeros(3, 1); mApplied = zeros(3, 1); Bold = [];
    for kk = 0:nSteps-1
        phase = mod(kk, N);
        if phase == mod(nOff, N)
            if nOff > 0 && tauS > 0, leakR = c*mHold*exp(-nOff/tauS);
            elseif nOff > 0, leakR = zeros(3, 1);
            else, leakR = c*mApplied; end
            Bmr = Btrue(kk + 1, :)' + leakR;
            if isempty(Bold), Bold = Btrue(1, :)'; end      % the model starts from the true field
            mHold = -k*(Bmr - Bold)/N;
            mHold = min(max(mHold, -P.mtqMax), P.mtqMax);
            Bold = Bmr;
        end
        if phase < nOff, mApplied = zeros(3, 1); else, mApplied = mHold; end
        mRef(kk + 1, :) = mApplied';
    end
    d = max(max(abs(mBlock - mRef)));
    fprintf('cycle %2d s, off %d s, c %.0e, tau %.1f s: max difference %.2e A m^2 (peak %.1f)\n', ...
        cases(i, 1), cases(i, 2), c, cases(i, 4), d, max(max(abs(mRef))));
    ok = ok && d < 1e-9;
end
if ok, fprintf('bdot_cycle_check: all cases agree\n'); else, fprintf('bdot_cycle_check: DIFFERENCE FOUND\n'); end
end
