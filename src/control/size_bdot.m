function s = size_bdot(p)
% Hand sizing of the B-dot gain for the command m = -k dB/dt.
% Averaged over an orbit that is near-polar to the dipole, |B|^2 averages 2.5 b0^2, and about
% 2/3 of it acts on a body turning about a random axis. The rate then decays as exp(-t/tauE).
a = p.Re + p.altitude;
period = 2*pi*sqrt(a^3/p.mu);
b0 = p.dipole/a^3;
jMean = trace(p.inertia)/3;
s.meanB2 = 2.5*b0^2;
s.tauE = p.bdotTauOrbits*period;
s.k = 1.5*jMean/(s.meanB2*s.tauE);
w0 = norm(p.tipOffRate);
s.tDetumble = s.tauE*log(w0/p.detumbleTarget);
s.mAtStart = s.k*w0*b0;
end
