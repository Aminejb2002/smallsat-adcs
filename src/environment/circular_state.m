function [r, v] = circular_state(a, inc, raan, u, p)
cO = cos(raan); sO = sin(raan);
ci = cos(inc); si = sin(inc);
cu = cos(u); su = sin(u);
r = a*[cu*cO - su*ci*sO; cu*sO + su*ci*cO; su*si];
v = sqrt(p.mu/a)*[-su*cO - cu*ci*sO; -su*sO + cu*ci*cO; cu*si];
end
