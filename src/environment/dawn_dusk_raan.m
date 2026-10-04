function raan = dawn_dusk_raan(daysJ2000)
s = sun_direction(daysJ2000);
raan = atan2(s(2), s(1)) - pi/2;
end
