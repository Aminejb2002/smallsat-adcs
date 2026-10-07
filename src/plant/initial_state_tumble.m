function ic = initial_state_tumble(p)
% Start after separation: same orbit as initial_state, arbitrary attitude, tip-off body rate.
base = initial_state(p, 0);
ic.r = base.r;
ic.v = base.v;
ic.q = p.initialQuat(:)/norm(p.initialQuat);
ic.w = p.tipOffRate(:);
end
