function wd = angular_accel(w, h, tauExt, tauW, inertia)
% Body angular acceleration with reaction wheels: w body rate, h wheel momentum, tauExt external
% torque on the body, tauW motor torque on the wheels (the body feels the reaction, -tauW).
w = w(:);
h = h(:);
wd = inertia \ (tauExt(:) - tauW(:) - cross(w, inertia*w + h));
end
