function tau = torque_srp(A, sunEci, lit, p)
force = -p.Psrp*p.srpCr*p.srpArea*(A*sunEci(:))*double(lit);
tau = cross(p.srpArm, force);
end
