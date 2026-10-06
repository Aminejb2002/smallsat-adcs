function [qHat, wHat, sig, bHat] = mekf_outputs(x, wMeas, qMeas)
% Estimates handed to the controller: attitude, rate, 1-sigma attitude uncertainty per axis (rad)
% and the gyro bias estimate. Before the first star fix the filter has not started: attitude and
% rate are the raw measurements, the bias is zero and sigma is zero.
if norm(x(1:4)) == 0
    qHat = qMeas(:);
    wHat = wMeas(:);
    sig = zeros(3, 1);
    bHat = zeros(3, 1);
    return
end
P = reshape(x(14:157), 12, 12);
qHat = x(1:4);
wHat = x(5:7);
sig = sqrt(diag(P(1:3, 1:3)));
bHat = x(8:10);
end
