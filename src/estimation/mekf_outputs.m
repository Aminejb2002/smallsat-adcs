function [qHat, wHat, sig, bHat] = mekf_outputs(x, wMeas, qMeas)
% Estimates handed to the controller: attitude, bias-corrected rate, 1-sigma attitude
% uncertainty per axis (rad) and the bias estimate. Before the first star fix the raw
% measurements are passed through.
if norm(x(1:4)) == 0
    qHat = qMeas(:);
    wHat = wMeas(:);
    sig = zeros(3, 1);
    bHat = zeros(3, 1);
    return
end
P = reshape(x(8:43), 6, 6);
qHat = x(1:4);
bHat = x(5:7);
wHat = wMeas(:) - bHat;
sig = sqrt(diag(P(1:3, 1:3)));
end
