function ic = initial_state_case(p, name)
% Initial states: 'tumble' after separation, 'point' after detumble with a large attitude error.
switch name
    case 'tumble'
        ic = initial_state_tumble(p);
    case 'point'
        ic = initial_state_tumble(p);
        ic.w = [0.05; -0.08; 0.06]*pi/180;
    otherwise
        error('initial_state_case:unknown', 'unknown case %s', name);
end
ic.h = zeros(3, 1);
end
