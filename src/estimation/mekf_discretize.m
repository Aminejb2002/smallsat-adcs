function [Phi, Qd] = mekf_discretize(F, Qc, dt)
% Transition matrix and process noise over dt, series expansions of the matrix exponential and
% of the Van Loan integral to second and third order in dt.
A = F*dt;
Phi = eye(size(F, 1)) + A + 0.5*(A*A);
FQ = F*Qc;
Qd = Qc*dt + (FQ + FQ')*dt^2/2 + (F*FQ + 2*(FQ*F') + (F*FQ)')*dt^3/6;
end
