function q = dcm_to_quat(A)
t = trace(A);
[~, k] = max([t, A(1,1), A(2,2), A(3,3)]);
switch k
    case 1
        q0 = 0.5*sqrt(1 + t);
        s = 4*q0;
        q = [q0; (A(2,3) - A(3,2))/s; (A(3,1) - A(1,3))/s; (A(1,2) - A(2,1))/s];
    case 2
        qx = 0.5*sqrt(1 + 2*A(1,1) - t);
        s = 4*qx;
        q = [(A(2,3) - A(3,2))/s; qx; (A(1,2) + A(2,1))/s; (A(1,3) + A(3,1))/s];
    case 3
        qy = 0.5*sqrt(1 + 2*A(2,2) - t);
        s = 4*qy;
        q = [(A(3,1) - A(1,3))/s; (A(1,2) + A(2,1))/s; qy; (A(2,3) + A(3,2))/s];
    otherwise
        qz = 0.5*sqrt(1 + 2*A(3,3) - t);
        s = 4*qz;
        q = [(A(1,2) - A(2,1))/s; (A(1,3) + A(3,1))/s; (A(2,3) + A(3,2))/s; qz];
end
if q(1) < 0
    q = -q;
end
end
