function Qc = mekf_process_noise(p, accelModel)
% Continuous process noise density of the error state [attitude; rate; bias; torque]. The rate
% gets white angular acceleration noise for unmodelled effects plus a term that grows with the
% modelled acceleration accelModel (rad/s^2, per axis): inertia and wheel torque are only known to
% p.mekfModelError, and that error acts over about p.mekfModelCorrTime seconds. The gyro bias and
% the disturbance torque are random walks.
qRate = p.mekfAccelNoise^2 + (p.mekfModelError*accelModel(:)).^2*p.mekfModelCorrTime;
Qc = blkdiag(zeros(3), diag(qRate), p.gyroRrw^2*eye(3), p.mekfTorqueRw^2*eye(3));
end
