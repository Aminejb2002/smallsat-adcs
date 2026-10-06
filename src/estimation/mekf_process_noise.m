function Qc = mekf_process_noise(p)
% Continuous process noise density of the error state: white angular acceleration on the rate
% (model error), gyro bias random walk and a random walk on the disturbance torque.
Qc = blkdiag(zeros(3), p.mekfAccelNoise^2*eye(3), p.gyroRrw^2*eye(3), p.mekfTorqueRw^2*eye(3));
end
