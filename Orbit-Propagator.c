#include <stdio.h>
#include <stdlib.h>
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
#define A1 3
#define A2 4
#define A3 5
#define mass 10

#define dist_to_sun 150e6

#define C_D 2.2

#define dt 1
#define total_time 100000

typedef struct{
    double x;
    double y;
    double z;
} vector;

typedef struct{
    vector x;
    vector y;
} accel;

double mag(vector v){
    return pow(pow(v.x,2) + pow(v.y,2) + pow(v.z,2),0.5);
}

vector scalar_mult(double scalar, vector v){
    vector mul_v = {scalar*v.x, scalar*v.y,scalar*v.z};
    return mul_v;
}

vector add(vector v1, vector v2){
    vector added = {v1.x+v2.x, v1.y+v2.y,v1.z+v2.z};
    return added;
}

vector a = {0, 0, 0};
vector velocity = {0, 7560, 0};
vector position = {0,3e6,7e6};

void print_vector(vector v){
    printf("%f,%f,%f\n", v.x,v.y,v.z);
}

vector get_gravity(double x, double y, double z){
    vector r = {x, y, z};
    double r_mag = mag(r);
    
    double ax = ((u*J2*(pow(Re, 2)))/2)*((15*x*pow(z, 2))/pow(r_mag, 7) - (3*x/pow(r_mag, 5)));
    double ay = ((u*J2*(pow(Re, 2)))/2)*((15*y*pow(z, 2))/pow(r_mag, 7) - (3*y/pow(r_mag, 5)));
    double az = ((u*J2*(pow(Re, 2)))/2)*((15*pow(z, 3))/pow(r_mag, 7) - (9*z/pow(r_mag, 5)));

    vector a_j2 = {ax, ay, az};
    vector a = scalar_mult(-((G*M_e)/(pow(r_mag, 3))), r);
    // return add(a_j2,a);
    return a;
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

double get_time(int julian_date){
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

    vector sat_sun_unit = scalar_mult(1/mag(sat_sun),sat_sun);
    vector r_unit = scalar_mult(1/mag(position),position);

    double theta = acos(dot(r_unit, earth_sun_unit));

    if ((theta > (pi/2)) && (mag(position)*sin(theta) < Re)) {
        vector zero = {0, 0, 0};
        return zero;
    }

    double P_weirdthing = S_weirdthing/(c*pow(mag(sat_sun), 2));

    double costheta1 = fabs(dot(normals[0], sat_sun_unit));
    double projArea1 = A1*costheta1;
    vector FSRP_face1 = scalar_mult(-P_weirdthing*projArea1, add(scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]), scalar_mult((1-R_spec), sat_sun_unit)));

//scalar_mult(-P_weirdthing*projArea1, add(scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]), scalar_mult((1-R_spec), sat_sun_unit)))
//add(scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]), scalar_mult((1-R_spec), sat_sun_unit))
//scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]);
//scalar_mult((1-R_spec), sat_sun_unit)


    double costheta2 = fabs(dot(normals[1], sat_sun_unit));
    double projArea2 = A2*costheta2;
    vector FSRP_face2 = scalar_mult(-P_weirdthing*projArea2, add(scalar_mult(2*(R_diff/3 + R_spec*costheta2), normals[1]), scalar_mult((1-R_spec), sat_sun_unit)));

    double costheta3 = fabs(dot(normals[2], sat_sun_unit));
    double projArea3 = A3*costheta3;
    vector FSRP_face3 = scalar_mult(-P_weirdthing*projArea3, add(scalar_mult(2*(R_diff/3 + R_spec*costheta3), normals[2]), scalar_mult((1-R_spec), sat_sun_unit)));

    vector FSRP = add(FSRP_face1, add(FSRP_face2, FSRP_face3));
    return FSRP;
}

vector drag(double rho, vector* normals, vector velocity){
    vector unit_velocity = scalar_mult(1/mag(velocity), velocity);
    // print_vector(unit_velocity);

    double costheta1 = fabs(dot(normals[0], unit_velocity));
    double projArea1 = A1*costheta1;
    double costheta2 = fabs(dot(normals[1], unit_velocity));
    double projArea2 = A2*costheta2;
    double costheta3 = fabs(dot(normals[2], unit_velocity));
    double projArea3 = A3*costheta3;
    if(projArea3 != projArea3){
        exit(0);
    }
    double A = projArea1 + projArea2 + projArea3;
    vector F_D = scalar_mult(((-rho*mag(velocity)*C_D*A)/2), velocity);
    // printf("rho %f\n", rho);
    return F_D;
}

double get_rho(vector position){
    int h0[] = {
        0, 25, 30, 35, 40, 45, 50, 55, 60, 65,
        70, 75, 80, 85, 90, 95, 100, 110, 120, 130,
        140, 150, 160, 180, 200, 250, 300, 350, 400, 450,
        500, 600, 700, 800, 900, 1000
    };

    double rho0[] = {
        1.225, 3.90E-02, 1.77E-02, 8.28E-03, 3.97E-03, 2.00E-03,
        1.06E-03, 5.82E-04, 3.21E-04, 1.72E-04, 8.77E-05, 4.18E-05,
        1.91E-05, 8.34E-06, 3.40E-06, 1.34E-06, 5.30E-07, 9.66E-08,
        2.44E-08, 8.48E-09, 3.85E-09, 2.07E-09, 1.22E-09, 5.46E-10,
        2.79E-10, 7.25E-11, 2.42E-11, 9.16E-12, 3.73E-12, 1.59E-12,
        6.97E-13, 1.45E-13, 3.61E-14, 1.17E-14, 5.25E-15, 3.02E-15
    };

    double H[] = {
        8.44, 6.49, 6.75, 7.07, 7.47, 7.83, 7.95, 7.73, 7.29, 6.81,
        6.33, 6.00, 5.70, 5.41, 5.38, 5.74, 6.15, 8.06, 11.60, 16.10,
        20.60, 24.60, 26.30, 33.20, 38.50, 46.90, 52.50, 56.40, 59.40,
        62.20, 65.80, 79.00, 109.00, 164.00, 225.00, 268.00
    };

    double h = mag(position) - Re;
    int i = 0;
    for (i = 0; i<36; i++) {
        if (h<=h0[i]*1000) break;
    }
    i--;

    if(i < 0){
        return 0;
    }

    double rho = rho0[i]*exp(-((h-(h0[i]*1000))/H[i]*1000));
    printf("rho %f\n",rho);
    printf("h %f\n",h);
    return rho;
}

int main(){
    FILE *fptr;
    fptr = fopen("Values.txt", "w");
    vector *pos = (vector*)malloc((total_time/dt)* sizeof(vector));
    for (int i = 1; i<=total_time; i+=dt){
        printf("mag %f\n", mag(position));
        vector a_t = get_gravity(position.x, position.y, position.z);
        vector *normals = (vector*)malloc(3* sizeof(vector));
        get_normals(normals);
        vector a_drag = drag(get_rho(position), normals, velocity);
        vector F_SRP = get_SRP(position, get_time(get_julian()), normals);
        vector a_SRP = scalar_mult(1/M_s, F_SRP);
        // a = add(add(a_t.x, a_t.y), add(a_drag, a_SRP));
        // a = add(a_t,a_drag);
        a = a_t;
        velocity = add(velocity, scalar_mult(dt, a));
        position = add(position, scalar_mult(dt, velocity));
        pos[(i-1)/dt] = position; 
        
        fprintf(fptr, "%f %f %f\n", position.x, position.y, position.z);
    }
    fclose(fptr);
}