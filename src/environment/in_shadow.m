function tf = in_shadow(r, s, Re)
% Cylindrical Earth shadow.
r = r(:);
along = r'*s(:);
tf = along < 0 && norm(r - along*s(:)) < Re;
end
