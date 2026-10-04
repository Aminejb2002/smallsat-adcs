function ref = reference_run(duration, step, tol, ic, ctrl)
% Plain-MATLAB propagation of orbit, attitude and reaction wheels, the golden check for the
% Simulink model. ctrl fields (all optional): gain and mMax for the sampled B-dot loop,
% pointing = true for the wheel pointing law (continuous, ideal state feedback).
% Noise-free magnetometer; the commanded dipole is held between sensor samples.
p = smallsat_params;
if nargin < 4 || isempty(ic)
    ic = initial_state(p);
end
if ~isfield(ic, 'h')
    ic.h = zeros(3, 1);
end
if nargin < 5 || isempty(ctrl)
    ctrl = struct();
end
if ~isfield(ctrl, 'gain'), ctrl.gain = 0; end
if ~isfield(ctrl, 'mMax'), ctrl.mMax = 0; end
if ~isfield(ctrl, 'pointing'), ctrl.pointing = false; end

Ts = p.sensorStep;
nStep = round(duration/Ts);
every = round(step/Ts);
if abs(nStep*Ts - duration) > 1e-9 || abs(every*Ts - step) > 1e-9
    error('reference_run:grid', 'duration and step must be multiples of the sensor step');
end

opts = odeset('RelTol', tol, 'AbsTol', tol);
x = [ic.r; ic.v; ic.q; ic.w; ic.h];
Bprev = body_field(0, x, p);

nRec = floor(nStep/every) + 1;
t = zeros(nRec, 1);
X = zeros(nRec, 16);
M = zeros(nRec, 3);
B = zeros(nRec, 3);
tau = zeros(nRec, 3);
parts = zeros(nRec, 12);
tauWheel = zeros(nRec, 3);
qRef = zeros(nRec, 4);
row = 0;
for k = 0:nStep
    tk = k*Ts;
    Bk = body_field(tk, x, p);
    m = -ctrl.gain*(Bk - Bprev)/Ts;
    m = min(max(m, -ctrl.mMax), ctrl.mMax);
    Bprev = Bk;
    if mod(k, every) == 0
        row = row + 1;
        [tauK, partsK] = disturbance_block(tk, x(1:3), x(4:6), x(7:10), p);
        t(row) = tk;
        X(row, :) = x';
        M(row, :) = m';
        B(row, :) = Bk';
        tau(row, :) = tauK';
        parts(row, :) = partsK(:)';
        tauWheel(row, :) = wheel_torque(x, p, ctrl.pointing)';
        qRef(row, :) = nadir_reference(x(1:3), x(4:6))';
    end
    if k < nStep
        [~, y] = ode113(@(tt, xx) rhs(tt, xx, p, m, ctrl.pointing), [tk, tk + Ts], x, opts);
        x = y(end, :)';
    end
end
ref.t = t;
ref.r = X(:, 1:3);
ref.v = X(:, 4:6);
ref.q = X(:, 7:10);
ref.w = X(:, 11:13);
ref.h = X(:, 14:16);
ref.m = M;
ref.B = B;
ref.tau = tau;
ref.parts = parts;
ref.tauW = tauWheel;
ref.qref = qRef;
end

function B = body_field(t, x, p)
B = quat_to_dcm(x(7:10))*magnetic_field_eci(x(1:3), t, p);
end

function tauW = wheel_torque(x, p, pointing)
if pointing
    [qref, wref] = nadir_reference(x(1:3), x(4:6));
    tauW = wheel_motor(-pointing_law(x(7:10), x(11:13), qref, wref, p), x(14:16), p);
else
    tauW = zeros(3, 1);
end
end

function dx = rhs(t, x, p, m, pointing)
tauDist = disturbance_block(t, x(1:3), x(4:6), x(7:10), p);
tauMt = cross(m, body_field(t, x, p));
tauW = wheel_torque(x, p, pointing);
wd = angular_accel(x(11:13), x(14:16), tauDist + tauMt, tauW, p.inertia);
dx = [x(4:6); j2_accel(x(1:3), p); quat_rate(x(7:10), x(11:13)); wd; tauW];
end
