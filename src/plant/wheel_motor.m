function tau = wheel_motor(cmd, h, p)
% Wheel motor torque: limited to the motor torque, and zero when it would push a wheel past its momentum limit.
tau = min(max(cmd(:), -p.wheelTorqueMax), p.wheelTorqueMax);
h = h(:);
blocked = (h >= p.wheelMomentumMax & tau > 0) | (h <= -p.wheelMomentumMax & tau < 0);
tau(blocked) = 0;
end
