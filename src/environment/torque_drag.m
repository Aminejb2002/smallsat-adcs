function tau = torque_drag(A, r, v, p)
r = r(:);
v = v(:);
vRel = v - cross([0; 0; p.wEarth], r);
h = norm(r) - p.Re;
force = -0.5*atmos_density(h, p)*p.dragCd*p.dragArea*norm(vRel)*(A*vRel);
tau = cross(p.dragArm, force);
end
