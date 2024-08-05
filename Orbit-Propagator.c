#include <stdio.h>
#include <math.h>

#define u 3.986e14
#define Re 6371837
#define J2  1.082629e-3
#define c  299792458
#define S_weirdthing  1362
#define R_diff  0.2
#define R_spec  0.5
#define R_abs  0.3
#define M_e  5.972e24
#define M_s 30
#define G 6.672e-11
#define pi 3.14

#define altitude 600000
#define A1  3
#define A2  4
#define A3 5

#define dist_to_sun 150e6

#define C_D = 2.2

#define dt 10

typedef struct{
    double x;
    double y;
    double z;
} vector;

double mag(vector v){
    return pow(pow(v.x,2) + pow(v.y,2) + pow(v.z,2),0.5);
}

vector scalar_mult(int scalar, vector v){
    vector mul_v = {scalar*v.x, scalar*v.y,scalar*v.z};
    return mul_v;
}

vector add(vector v1, vector v2){
    vector added = {v1.x+v2.x, v1.y+v2.y,v1.z+v2.z};
    return added;
}

vector get_gravity(vector position, vector velocity){
    double r_mag = mag(position);
    double ax = u * J2 * (pow(Re,2)/2) * ((15*position.x*pow(position.z,2))/pow(r_mag,7) - (3*position.x/pow(r_mag,5)));     
    double ay = u * J2 * (pow(Re,2)/2) * ((15*position.y*pow(position.z,2))/pow(r_mag,7) - (3*position.y/pow(r_mag,5)));     
    double az = u * J2 * (pow(Re,2)/2) * ((15*pow(position.z,3))/pow(r_mag,7) - (9*position.z/pow(r_mag,5)));

    vector a = {ax,ay,az};
    vector g = scalar_mult(G*M_e/pow(r_mag,3),position);

    return add(a,g);
}

void get_normals(vector* normals){
    vector n1 = {1,0,0};
    vector n2 = {0,1,0};
    vector n3 = {0,0,1};

    normals[0] = n1;
    normals[1] = n2;
    normals[2] = n3;
}

int get_julian(){
    return 2600000;
}

double get_time(julian_date){
    return (julian_date - 2451545)/36525;
}

double dot(vector v1, vector v2){
    return v1.x*v2.x + v1.y*v2.y + v1.z* v2.z;
}

vector get_SRP(vector position,double time,vector* normals){
    double phi = 280.460 + 36000.771*time;
    double M = 357.5277233 + 35999.05034*time; 
    double phi_ecl = phi + 1.914666471*sin(M) + 0.019994643*sin(2*M); 
    double epsilon = 23.439291 - 0.0130042*time;

    vector earth_sun_unit = {cos(phi_ecl), cos(epsilon)*sin(phi_ecl), sin(epsilon)*sin(phi_ecl)};
    vector sat_sun = add(scalar_mult(dist_to_sun*epsilon,earth_sun_unit),scalar_mult(-1,position));

    vector sat_sun_unit = scalar_mult(1/norm(sat_sun),sat_sun);
    vector r_unit = scalar_mult(1/norm(position),position);

    double theta = acos(dot(r_unit, earth_sun_unit));

    if ((theta > (pi/2)) && (norm(position)*sin(theta) < Re)) {
        vector zero = {0, 0, 0};
        return zero;
    }

    double P_weirdthing = S_weirdthing/(c*pow(norm(sat_sun), 2));

    double costheta1 = abs(dot(normals[0], sat_sun_unit));
    double projArea1 = A1*costheta1;
    vector FSRP_face1 = -P_weirdthing*projArea1


    double costheta2 = abs(dot(normals[1], sat_sun_unit));
    double projArea2 = A2*costheta2;

    double costheta3 = abs(dot(normals[2], sat_sun_unit));
    double projArea3 = A3*costheta3;
}