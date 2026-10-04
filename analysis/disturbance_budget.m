function b = disturbance_budget
% Disturbance torques over one orbit in nadir pointing, with hand-estimate bounds.

p = smallsat_params;
a = p.Re + p.altitude;
n = sqrt(p.mu/a^3);
period = 2*pi/n;
inc = sso_inclination(a, p);
raan = dawn_dusk_raan(p.epochDays);
sun = sun_direction(p.epochDays);
hHat = [sin(inc)*sin(raan); -sin(inc)*cos(raan); cos(inc)];

t = 0:1:period;
names = {'gravity gradient', 'solar pressure', 'drag', 'magnetic'};
tau = zeros(numel(t), 3, 4);
shadowSamples = 0;
for i = 1:numel(t)
    [r, v] = circular_state(a, inc, raan, n*t(i), p);
    A = nadir_dcm(r, v);
    lit = ~in_shadow(r, sun, p.Re);
    shadowSamples = shadowSamples + ~lit;
    tau(i, :, 1) = torque_gravity_gradient(A*r, p.inertia, a, p)';
    tau(i, :, 2) = torque_srp(A, sun, lit, p)';
    tau(i, :, 3) = torque_drag(A, r, v, p)';
    tau(i, :, 4) = cross(p.residualDipole, A*magnetic_field_eci(r, t(i), p))';
end

b0 = p.dipole/a^3;
vOrb = sqrt(p.mu/a);
spread = max(diag(p.inertia)) - min(diag(p.inertia));
hand = [1.5*p.mu/a^3*spread, ...
        p.Psrp*p.srpCr*p.srpArea*norm(p.srpArm), ...
        0.5*p.rhoRef*p.dragCd*p.dragArea*vOrb^2*norm(p.dragArm), ...
        norm(p.residualDipole)*2*b0];

stats = repmat(torque_stats(tau(:, :, 1), period), 1, 5);
for k = 1:4
    stats(k) = torque_stats(tau(:, :, k), period);
end
stats(5) = torque_stats(sum(tau, 3), period);

b.names = [names, {'total'}];
b.periodSeconds = period;
b.inclinationDeg = inc*180/pi;
b.betaDeg = asind(hHat'*sun);
b.shadowSamples = shadowSamples;
b.samples = numel(t);
b.peak = [stats.peak];
b.hand = [hand, sum(hand)];
b.meanNorm = [stats.meanNorm];
b.secularPerOrbit = [stats.secular];
b.cyclicSwing = [stats.swing];
b.secularPerDay = stats(5).meanNorm*86400;

if nargout == 0
    fprintf('altitude %.0f km, period %.1f s, SSO inclination %.3f deg\n', ...
        p.altitude/1e3, b.periodSeconds, b.inclinationDeg);
    fprintf('beta %.1f deg, shadow samples %d of %d\n', b.betaDeg, b.shadowSamples, b.samples);
    fprintf('%-18s %17s %12s %12s %16s %15s\n', 'source', 'peak |tau| [N m]', 'hand bound', ...
        'mean |tau|', 'secular H/orbit', 'cyclic swing H');
    for k = 1:5
        fprintf('%-18s %17.3e %12.3e %12.3e %16.3e %15.3e\n', b.names{k}, b.peak(k), b.hand(k), ...
            b.meanNorm(k), b.secularPerOrbit(k), b.cyclicSwing(k));
    end
    fprintf('secular H per day [N m s]: %.4f\n', b.secularPerDay);
    clear b
end
end

function s = torque_stats(tt, period)
m = mean(tt, 1);
cum = cumsum(tt - m, 1);
s.peak = max(vecnorm(tt, 2, 2));
s.meanNorm = norm(m);
s.secular = norm(m)*period;
s.swing = max(max(cum, [], 1) - min(cum, [], 1))/2;
end
