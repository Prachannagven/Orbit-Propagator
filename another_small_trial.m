% Sample table (altitude in meters, orbital velocity in m/s)
altitude_table = [100e6; 200e6; 300e6; 400e6; 500e6];
orbital_velocity_table = [7500; 7600; 7700; 7800; 7900];

% Lookup altitude (replace with your desired value)
lookup_altitude = 150;

% Perform linear interpolation
lookup_velocity = interp1(array_h0, array_rho0, lookup_altitude);

disp(['Orbital velocity at ', num2str(lookup_altitude), ' meters: ', num2str(lookup_velocity), ' m/s']);
