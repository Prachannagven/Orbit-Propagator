%Constants
constants.u = 3.986 * 10^14;
constants.J2 = 1.082629*10^-2;
constants.Re = 6378137;
constants.c = 299792458;
constants.S_weirdthing = 1362;
constants.R_diff = 0.2;
constants.R_spec = 0.5;
constants.R_abs = 0.3;

constants.altitude = 600000;
constants.rho_0 = []; 

constants.A1 = 3;
constants.A2 = 4;
constants.A3 = 5;

constants.dist_to_sun = 150*10^6;

constants.C_D = 2.2; %default value of GMAT model

%Why Bro, Why

a = [0,0,0];
v = [0,0,0];
x = 500000;
y = 400000;
z = 50000;
r = [x,y,z];
mass = 10;

dt = 0.01;

x = getNormals();
for i = 1:1:1000
    a_g = acc(x,y,z,constants)
    
end

%Functions 

%For x,y,z acceleration
function [A] = projArea(n1, n2, n3, v) %to get the area along the direction of motion
    A = A1*(abs(dot(n1, v))) + A2*(abs(dot(n2, v))) + A3*(abs(dot(n3, v)));
end

function a = acc(x,y,z,constants) %to get the gravitational acceleration of the body over time
    r = sqrt(x^2+y^2+z^2);
    ax = ((constants.u*constants.J2*((constants.Re)^2))/2)*((15*x*z^2)/r^7 - (3*x/r^5));
    ay = ((constants.u*constants.J2*(constants.Re^2))/2)*((15*y*z^2)/r^7 - (3*y/r^5));
    az = ((constants.u*constants.J2*(constants.Re^2))/2)*((15*z^3)/r^7 - (9*z/r^5));

    a = [ax,ay,az];
end

function [rho] = rho(h) %to get the value of atmospheric density at differnet time
    %add lookup functions for h_0, H and rho_0 and then solve
    rho = rho_0*exp(-((h-h_0)/H));
end

function [n1, n2, n3] = getNormals() %getting the nromals from the whtaevers as input (find out what whatever is)
    n1(3) = [1, 0, 0];
    n2(3) = [0, 1, 0];
    n3(3) = [0, 0, 1];
end

function [F_D] = drag(rho, n1, n2, n3, v) %getting drag force
    A = projArea(n1, n2, n3, v);
    F_D = -((rho*v^2*C_D*A)/2)*(v ./ norm(v));
end

function [JulianDate] = getJulian()
    JulianDate = 2600000;
end

function [T] = getTime(JulianDate)
    T = (JulianDate - 2451545)/36525;
end


function [F_SRP] = solarPressure(x, y, z, T, n1, n2, n3)
    r(3) = [x, y z]; %a prereq variable to simplify stuff

    %calculating the constants we need to even actually get solar pressure
    phi = 280.460 + 36000.771*T; %mean longitude of sun
    M = 357.5277233 + 35999.05034*T; %mean anomaly of sun
    phi_ecl = phi + 1.914666471*sin(M) + 0.019994643*sin(2*M); %longtitude of the cliptic
    epsilon = 23.439291 - 0.0130042*T; %obliquity of the ecliptic

    %getting earth to sun vector
    earth_sun_unit(3) = [cos(phi_ecl), cos(epsilon)*sin(phi_ecl), sin(epislon)*sin(phi_ecl)];
    sat_sun = dist_to_sun*epsilon_ES - r;
    sat_sun_unit = sat_sun ./ norm(sat_sun);
    r_unit = r ./ norm(r);
    theta = acos(dot(r_unit, earth_sun_unit));

    if(theta > (pi/2) && (norm(r)*sin(theta)<Re))
        return;
    end

    %One more step before doing it
    P_weirdthing = S_weirdthing/(c*(norm(sat_sun))^2);

    %now to calculate piece by piece
    costheta1 = abs(dot(n1, sat_sun_unit));
    projArea1 = A1*costheta1;
    FSRP_face1 = -P_weirdthing*projArea1*(2*n1*(R_diff/3 + R_spec*costheta1) + (1-R_spec)*sat_sun_unit);

    costheta2 = abs(dot(n2, sat_sun_unit));
    projArea2 = A2*costheta2;
    FSRP_face2 = -P_weirdthing*projArea2*(2*n2*(R_diff/3 + R_spec*costheta2) + (1-R_spec)*sat_sun_unit);

    costheta3 = abs(dot(n3, sat_sun_unit));
    projArea3 = A3*costheta3;
    FSRP_face3 = -P_weirdthing*projArea3*(2*n3*(R_diff/3 + R_spec*costheta3) + (1-R_spec)*sat_sun_unit);

    FSRP = FSRP_face1 + FSRP_face2 +FSRP_face3;
end


