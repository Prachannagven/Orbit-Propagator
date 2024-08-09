function satellite_propagation(initial_altitude, mass, dimensions)

% Constants
G = 6.67430e-11; % Gravitational constant (m^3 kg^-1 s^-2)
M_earth = 5.972e24; % Mass of Earth (kg)
R_earth = 6371000; % Earth's radius (m)
J2 = 1.08263e-3; % Earth's J2 coefficient
Cd = 2.2; % Drag coefficient (assumed)
A = 0.01 * 0.03; % Reference area for drag (m^2)
rho0 = 1.225; % Air density at sea level (kg/m^2)
H = 7200; % Scale height of the atmosphere (m)
P_solar = 1361; % Solar constant (W/m^2)
C_r = 1.2; % Solar radiation pressure coefficient

% Initial conditions
altitude = initial_altitude;
r = R_earth + altitude;
v = sqrt(G * M_earth / r); % Initial circular velocity

% Satellite state vector (position and velocity)
state = [r * [cos(0), sin(0), 0]; v * [-sin(0), cos(0), 0]];

% Simulation parameters
dt = 60; % Time step (seconds)
simulation_time = 86400; % Simulation time (seconds, one day)

% Time vector
t = 0:dt:simulation_time;

% Pre-allocate state matrix
states = zeros(length(t), 6);
states(1,:) = state;

for i = 2:length(t)
  % Calculate position and velocity
  r = state(i-1, 1:3);
  v = state(i-1, 4:6);

  % Calculate accelerations
  acc_gravity = -G * M_earth / norm(r)^3 * r;
  acc_j2 = j2_perturbation(r, J2, M_earth, R_earth);
  acc_drag = drag_acceleration(r, v, Cd, A, rho0, H);
  acc_srp = srp_acceleration(r, P_solar, A, C_r);
  acc_total = acc_gravity + acc_j2 + acc_drag + acc_srp;

  % Update state using numerical integration (e.g., Euler's method)
  state(i, 1:3) = state(i-1, 1:3) + dt * state(i-1, 4:6);
  state(i, 4:6) = state(i-1, 4:6) + dt * acc_total;
end

% Plot results
plot(states(:,1), states(:,2));
xlabel('X');
ylabel('Y');
title('Satellite Orbit');
end

function acc_j2 = j2_perturbation(r, J2, M_earth, R_earth)
  % J2 perturbation calculation (as before)
end

function acc_drag = drag_acceleration(r, v, Cd, A, rho0, H)
  % Calculate air density based on altitude
  altitude = norm(r) - R_earth;
  rho = rho0 * exp(-altitude / H);

  % Calculate drag force
  drag_force = -0.5 * rho * norm(v)^2 * Cd * A * v / norm(v);

  % Calculate drag acceleration
  acc_drag = drag_force / mass;
end

function acc_srp = srp_acceleration(r, P_solar, A, C_r)
  % Calculate solar radiation pressure acceleration
  AU = 149.6e9; % Astronomical unit (m)
  solar_flux = P_solar / (4 * pi * AU^2);
  srp_pressure = solar_flux * C_r;
  acc_srp = -srp_pressure * A / mass * r / norm(r);
end
