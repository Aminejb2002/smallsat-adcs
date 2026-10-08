function sc = globe_intro()
% Opens the 3D globe with the orbit of the simulated satellite (600 km, Sun-synchronous, dawn-dusk) so it can be
% screen-recorded as an intro clip. This is a plain two-body orbit for the picture only; the simulation itself
% uses the J2 model in src/environment. Needs Aerospace Toolbox and Satellite Communications Toolbox.
re = 6378137;
a = re + 600e3;
inc = 97.788;   % Sun-synchronous at 600 km, see tests/reference_check
start = datetime(2026, 6, 21, 12, 0, 0, 'TimeZone', 'UTC');
period = 2*pi*sqrt(a^3/3.986004418e14);
% right ascension of the ascending node for a dawn-dusk plane: 90 deg from the Sun's right ascension
n = days(start - datetime(2000, 1, 1, 12, 0, 0, 'TimeZone', 'UTC'));
L = 280.460 + 0.9856474*n;
g = 357.528 + 0.9856003*n;
lam = L + 1.915*sind(g) + 0.020*sind(2*g);
eps = 23.439;
sunRA = atan2d(cosd(eps)*sind(lam), cosd(lam));
raan = mod(sunRA + 90, 360);
sc = satelliteScenario(start, start + seconds(2*period), 10);
sat = satellite(sc, a, 0, inc, raan, 0, 0, 'Name', 'smallsat');
groundTrack(sat, 'LeadTime', period/2, 'TrailTime', period/2);
v = satelliteScenarioViewer(sc, 'CameraReferenceFrame', 'Inertial');
fprintf('Scenario ready. Orbit period %.0f s, RAAN %.1f deg.\n', period, raan);
fprintf('Start your screen recorder, then run: play(sc, ''PlaybackSpeedMultiplier'', 600)\n');
end
