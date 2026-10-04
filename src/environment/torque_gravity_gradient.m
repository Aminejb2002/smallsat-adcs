function tau = torque_gravity_gradient(rBody, inertia, rn, p)
rh = rBody(:)/rn;
tau = 3*p.mu/rn^3*cross(rh, inertia*rh);
end
