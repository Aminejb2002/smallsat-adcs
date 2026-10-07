function dump_check
% Checks of the magnetic momentum dumping law dump_law.
p = smallsat_params;
rng(5);
fails = 0;
total = 0;
worstPerp = 0; worstTorque = 0; worstSign = -Inf; worstScale = 0;
for k = 1:200
    B = (0.2 + 4*rand)*3e-5*randn(3, 1)/1;
    h = 0.05*randn(3, 1);
    p.dumpGain = 1/900;
    m = dump_law(B, h, p);
    tau = cross(m, B);
    want = -p.dumpGain*h;
    wantPerp = want - B*(B'*want)/(B'*B);
    unsat = max(abs(cross(B, want)/(B'*B))) <= p.mtqMax;
    if unsat
        worstTorque = max(worstTorque, norm(tau - wantPerp)/max(norm(wantPerp), eps));
    end
    worstPerp = max(worstPerp, abs(m'*B)/(norm(m)*norm(B) + eps));
    worstSign = max(worstSign, tau'*h/(norm(tau)*norm(h) + eps));
    worstScale = max(worstScale, max(abs(m))/p.mtqMax);
end
[fails, total] = report(fails, total, 'torque = part of -kh*h perpendicular to B', worstTorque, 1e-9, worstTorque < 1e-9);
[fails, total] = report(fails, total, 'dipole perpendicular to field', worstPerp, 1e-9, worstPerp < 1e-9);
[fails, total] = report(fails, total, 'torque never raises |h| (cos angle to h)', worstSign, 0, worstSign <= 1e-12);
[fails, total] = report(fails, total, 'dipole within the torquer limit', worstScale - 1, 1e-12, worstScale <= 1 + 1e-12);
m0 = dump_law([1e-5; -2e-5; 3e-5], zeros(3, 1), p);
[fails, total] = report(fails, total, 'zero momentum gives zero dipole', norm(m0), 1e-15, norm(m0) < 1e-15);
big = dump_law([3e-5; 0; 0], [0; 5; 5], p);
d1 = big/norm(big);
ref = cross([3e-5; 0; 0], [0; -5; -5]);
ref = ref/norm(ref);
[fails, total] = report(fails, total, 'saturation keeps the direction', norm(d1 - ref), 1e-12, norm(d1 - ref) < 1e-12);
fprintf('\n%d of %d checks passed\n', total - fails, total);
if fails > 0
    error('dump_check:failed', '%d checks failed', fails);
end
end

function [fails, total] = report(fails, total, name, value, limit, ok)
total = total + 1;
if ok, s = 'PASS'; else, s = 'FAIL'; fails = fails + 1; end
fprintf('%s  %-44s %10.3e (limit %.0e)\n', s, name, value, limit);
end
