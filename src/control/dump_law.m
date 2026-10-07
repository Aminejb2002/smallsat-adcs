function m = dump_law(B, h, p)
% Magnetic wheel momentum dumping. The wheels keep the attitude, so an external torque tau on the
% body changes the wheel momentum by the same amount; tau = -kh*h therefore drives h to zero.
% The torquer can only produce the part of tau perpendicular to B: m = B x tau / |B|^2.
% m is scaled (direction kept) so no component exceeds p.mtqMax.
B = B(:);
tau = -p.dumpGain*h(:);
m = cross(B, tau)/max(B'*B, eps);
peak = max(abs(m));
if peak > p.mtqMax
    m = m*(p.mtqMax/peak);
end
end
