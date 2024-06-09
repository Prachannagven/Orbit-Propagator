clc
clear

%Constants
constants.u = 3.986 * 10^14;
constants.J2 = 1.082629*10^-3;
constants.Re = 6371837;
constants.c = 299792458;
constants.S_weirdthing = 1362;
constants.R_diff = 0.2;
constants.R_spec = 0.5;
constants.R_abs = 0.3;
constants.M_e = 5.972e24;
constants.M_s = 30;
constants.G = 6.672e-11;

constants.altitude = 600000;
constants.rho_0 = []; 

constants.A1 = 3;
constants.A2 = 4;
constants.A3 = 5;

constants.dist_to_sun = 150*10^6;

constants.C_D = 2.2; %default value of GMAT model

%extracting table data so that we can use it in the drag force
data_table = readtable("Drag_Constants.xlsx", Sheet="Sheet1");

%Why Bro, Why

a = [0;0;0];
v = [0;7560;0];
x = 7e6;
y = 0;
z = 0;
r = [x;y;z];
mass = 10;

pos = [r];

dt = 10;
total_time = 1e5;


overall = [];

for j = 1:dt:total_time
    [a_g, a_J2] = acc(r(1),r(2),r(3),constants);
    a_drag = drag(rho(norm(r), data_table), [[1;0;0], [0;1;0], [0;0;1]], v, constants);
    F_SRP = solarPressure(r(1), r(2), r(3), getTime(getJulian()), [[1;0;0], [0;1;0], [0;0;1]], constants);
    a_SRP = F_SRP / constants.M_s;
    a = a_g + a_J2 + a_drag + a_SRP;
    v = v + a*dt;
    r = r + v*dt;
    pos = [pos r];
end

comet3(pos(1,:), pos(2,:), pos(3,:));

%Functions 

%For x,y,z acceleration
function [A] = projArea(normals, v, constants) %to get the area along the direction of motion
    A = constants.A1*(abs(dot(normals(1:3,1), v))) + constants.A2*(abs(dot(normals(1:3,2), v))) + constants.A3*(abs(dot(normals(1:3,3), v)));
end

function [a, a_j2] = acc(x,y,z,constants) %to get the gravitational acceleration of the body over time
    r = [x;y;z];
    r_mag = norm(r);
    ax = ((constants.u*constants.J2*((constants.Re)^2))/2)*((15*x*z^2)/r_mag^7 - (3*x/r_mag^5));
    ay = ((constants.u*constants.J2*(constants.Re^2))/2)*((15*y*z^2)/r_mag^7 - (3*y/r_mag^5));
    az = ((constants.u*constants.J2*(constants.Re^2))/2)*((15*z^3)/r_mag^7 - (9*z/r_mag^5));
    
    a_j2 = [ax;ay;az];
    a = -(constants.G*constants.M_e/(r_mag^3))*r;
end

function [rho] = rho(h_0, table_pull) %to get the value of atmospheric density at differnet time
    %getting individual double arrays so that we can use the interpolate function
    array_h0 = table2array([table_pull(:,1)]);
    
    array_rho0 = table2array([table_pull(:,2)]);
    array_H = table2array([table_pull(:,3)]);

    %using the interp1 function in order to get the values of rho0 and H
    H = interp1(array_h0, array_H, h_0);
    rho_0 = interp1(array_h0, array_rho0, h_0);

    %plugging in the values of rho0 and H so that we can get the value of rho
    rho = rho_0*exp(-((h-h_0)/H));
end

function [n1,n2,n3] = getNormals() %getting the nromals from the whtaevers as input (find out what whatever is)
    %for now just setting normals as constants and assuming the bodydoesn't rotate
    n1 = [1; 0; 0];
    n2 = [0; 1; 0];
    n3 = [0; 0; 1];
end

function [F_D] = drag(rho, normals, v, constants) %getting drag force
    A = projArea([normals(1:3,1), normals(1:3,2), normals(1:3,3)], v, constants);
    F_D = -((rho*(norm(v))^2*constants.C_D*A)/2)*(v ./ norm(v));
end

function [JulianDate] = getJulian()
    JulianDate = 2600000;
end

function [T] = getTime(JulianDate)
    T = (JulianDate - 2451545)/36525;
end


function [FSRP] = solarPressure(x, y, z, T, normals, constants)
    
    r = [x; y; z]; %a prereq variable to simplify stuff
    FSRP = 0;

    %calculating the constants we need to even actually get solar pressure
    phi = 280.460 + 36000.771*T; %mean longitude of sun
    M = 357.5277233 + 35999.05034*T; %mean anomaly of sun
    phi_ecl = phi + 1.914666471*sin(M) + 0.019994643*sin(2*M); %longtitude of the cliptic
    epsilon = 23.439291 - 0.0130042*T; %obliquity of the ecliptic

    %getting earth to sun vector
    earth_sun_unit = [cos(phi_ecl); cos(epsilon)*sin(phi_ecl); sin(epsilon)*sin(phi_ecl)];

    %What is epsilon_ES
    sat_sun = constants.dist_to_sun*epsilon*earth_sun_unit - r;

    sat_sun_unit = sat_sun ./ norm(sat_sun);
    r_unit = r ./ norm(r);
    theta = acos(dot(r_unit, earth_sun_unit));

    if(theta > (pi/2) && (norm(r)*sin(theta)<constants.Re))
        return;
    end

    %One more step before doing it
    P_weirdthing = constants.S_weirdthing/(constants.c*(norm(sat_sun))^2);

    %now to calculate piece by piece
    costheta1 = abs(dot(normals(1:3,1), sat_sun_unit));
    projArea1 = constants.A1*costheta1;
    FSRP_face1 = -P_weirdthing*projArea1*(2*normals(1)*(constants.R_diff/3 + constants.R_spec*costheta1) + (1-constants.R_spec)*sat_sun_unit);

    costheta2 = abs(dot(normals(1:3,2), sat_sun_unit));
    projArea2 = constants.A2*costheta2;
    FSRP_face2 = -P_weirdthing*projArea2*(2*normals(2)*(constants.R_diff/3 + constants.R_spec*costheta2) + (1-constants.R_spec)*sat_sun_unit);

    costheta3 = abs(dot(normals(1:3,3), sat_sun_unit));
    projArea3 = constants.A3*costheta3;
    FSRP_face3 = -P_weirdthing*projArea3*(2*normals(3)*(constants.R_diff/3 + constants.R_spec*costheta3) + (1-constants.R_spec)*sat_sun_unit);

    FSRP = FSRP_face1 + FSRP_face2 +FSRP_face3;
end


