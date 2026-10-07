function res = run_dumping(orbits, gains)
% Wheel momentum build-up under the pointing law with and without magnetic dumping. Plain
% propagation (fixed-step RK4, 1 s) with ideal state feedback, noise-free magnetometer: this
% isolates the dumping law from the estimator. The dipole is held between 1 s samples.
% gains: dump gains in 1/s, 0 = no dumping. Reports the wheel momentum, the magnetic torque and
% the pointing error.
if nargin < 1, orbits = 5; end
if nargin < 2, gains = [0 1/3600 1/1800 1/900]; end
p0 = smallsat_params;
a = p0.Re + p0.altitude;
period = 2*pi*sqrt(a^3/p0.mu);
T = orbits*period;
res = struct('gain', num2cell(gains));
fprintf('%9s %10s %10s %10s %10s %10s %10s %10s\n', 'dump tau', 'h final', 'h peak', 'h mean', 'tau_m rms', 'm peak', 'sat', 'point err');
fprintf('%9s %10s %10s %10s %10s %10s %10s %10s\n', '[s]', '[Nms]', '[Nms]', '[Nms]', '[Nm]', '[Am2]', 'fraction', 'mean [deg]');
for c = 1:numel(gains)
    p = p0;
    p.dumpGain = gains(c);
    ic = initial_state(p);
    x = [ic.r; ic.v; ic.q; ic.w; zeros(3, 1)];
    nStep = round(T);
    hN = zeros(nStep, 1); tm = zeros(nStep, 1); mP = zeros(nStep, 1); err = zeros(nStep, 1);
    for k = 0:nStep - 1
        t = k;
        B = quat_to_dcm(x(7:10))*magnetic_field_eci(x(1:3), t, p);
        m = zeros(3, 1);
        if p.dumpGain > 0
            m = dump_law(B, x(14:16), p);
        end
        f = @(tt, xx) rhs(tt, xx, p, m);
        k1 = f(t, x); k2 = f(t + 0.5, x + 0.5*k1); k3 = f(t + 0.5, x + 0.5*k2); k4 = f(t + 1, x + k3);
        x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
        x(7:10) = x(7:10)/norm(x(7:10));
        hN(k + 1) = norm(x(14:16));
        tm(k + 1) = norm(cross(m, B));
        mP(k + 1) = max(abs(m));
        [qref, ~] = nadir_reference(x(1:3), x(4:6));
        qe = quat_mul([qref(1); -qref(2:4)], x(7:10));
        err(k + 1) = 2*atan2(norm(qe(2:4)), abs(qe(1)))*180/pi;
    end
    sel = round(nStep/5) + 1:nStep;
    res(c).hFinal = hN(end); res(c).hPeak = max(hN); res(c).hMean = mean(hN(sel));
    res(c).tauM = sqrt(mean(tm.^2)); res(c).mPeak = max(mP);
    res(c).sat = mean(mP >= p.mtqMax*(1 - 1e-9)); res(c).err = mean(err(sel));
    if gains(c) > 0, tauS = 1/gains(c); else, tauS = Inf; end
    fprintf('%9.0f %10.4f %10.4f %10.4f %10.2e %10.2f %10.3f %10.4f\n', tauS, hN(end), max(hN), mean(hN(sel)), res(c).tauM, res(c).mPeak, res(c).sat, res(c).err);
end
end

function dx = rhs(t, x, p, m)
tauDist = disturbance_block(t, x(1:3), x(4:6), x(7:10), p);
B = quat_to_dcm(x(7:10))*magnetic_field_eci(x(1:3), t, p);
tauMt = cross(m, B);
[qref, wref] = nadir_reference(x(1:3), x(4:6));
w = x(11:13); h = x(14:16);
cmd = -pointing_law(x(7:10), w, qref, wref, p) - cross(w, p.inertia*w + h);
tauW = wheel_motor(cmd, h, p);
wd = angular_accel(w, h, tauDist + tauMt, tauW, p.inertia);
dx = [x(4:6); j2_accel(x(1:3), p); quat_rate(x(7:10), w); wd; tauW];
end
