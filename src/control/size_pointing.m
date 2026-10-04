function s = size_pointing(p)
% Gains of the pointing law and hand estimates of the slew from a 180 degree error.
% Per axis, kp = J wn^2 and kd = 2 zeta wn J (a second-order loop with a rate-limited approach).
j = diag(p.inertia);
s.wn = p.pointWn;
s.zeta = p.pointZeta;
s.wSlew = p.slewRateMax;
s.kp = j*s.wn^2;
s.kd = 2*s.zeta*s.wn*j;
s.thetaLinear = 2*s.zeta*s.wSlew/s.wn;
s.tolerance = 0.1*pi/180;
s.tCapture = (pi - s.thetaLinear)/s.wSlew + log(s.thetaLinear/s.tolerance)/(s.zeta*s.wn);
s.hSlew = max(j)*s.wSlew;
end
