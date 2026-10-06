function R = star_covariance(p)
% Star tracker attitude noise covariance in the body frame, rad^2: starSigma across the boresight
% and starBoresightRatio*starSigma about it.
b = p.starBoresight(:)/norm(p.starBoresight);
R = p.starSigma^2*(eye(3) + (p.starBoresightRatio^2 - 1)*(b*b'));
end
