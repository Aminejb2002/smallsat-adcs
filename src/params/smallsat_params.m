function p = smallsat_params
% Parameters of the generic Earth-observation smallsat and the environment.
% Spacecraft values are engineering assumptions, not data of a real satellite.

p.mu = 3.986004418e14;
p.Re = 6378.137e3;
p.J2 = 1.08262668e-3;
p.wEarth = 7.2921159e-5;
p.yearSeconds = 365.2422*86400;
p.Psrp = 4.56e-6;
p.dipole = 7.94e15;
p.epochDays = days(datetime(2026, 3, 20, 12, 0, 0) - datetime(2000, 1, 1, 12, 0, 0));

p.mass = 150;
p.inertia = [26.0 0.6 -0.4; 0.6 24.0 0.5; -0.4 0.5 16.0];
p.altitude = 600e3;

p.srpArea = 2.4;
p.srpCr = 1.3;
p.srpArm = [0.05; 0.02; 0.03];
p.dragArea = 1.5;
p.dragCd = 2.2;
p.dragArm = [0.03; -0.02; 0.04];
p.rhoRef = 1.0e-13;
p.scaleHeight = 70e3;
p.residualDipole = [0.3; -0.2; 0.25];

p.sensorStep = 1;
p.magSigma = 1e-7;
p.mtqMax = 30;
p.tipOffRate = [3; -2; 4]*pi/180;
p.detumbleTarget = 0.5*pi/180;
p.detumbleExit = 2.0*pi/180;
p.bdotTauOrbits = 0.5;

p.wheelTorqueMax = 0.05;
p.wheelMomentumMax = 2.0;
p.pointWn = 0.05;
p.pointZeta = 0.9;
p.slewRateMax = 1.0*pi/180;

p.gyroStep = 0.1;
p.gyroArw = 0.003*pi/180;
p.gyroRrw = 1.0e-9;
p.gyroBias0 = [1.0; -0.7; 0.5]*pi/180/3600;
p.starStep = 1;
p.starSigma = 10*pi/180/3600;
p.starRateLimit = 2.0*pi/180;

p.mekfBiasSigma0 = 1.0e-5;
end
