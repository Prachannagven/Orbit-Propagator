#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define u 3.986e14
#define R_e 6371837
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
#define mass 10.0

#define dist_to_sun 150e6

#define C_D 2.2

#define dt 10
#define total_time 750000

typedef struct{
    long double x;
    long double y;
    long double z;
} vector;

long double mag(vector v){
    return pow(pow(v.x,2) + pow(v.y,2) + pow(v.z,2),0.5);
}

vector scalar_mult(long double scalar, vector v){
    vector mul_v = {scalar*v.x, scalar*v.y,scalar*v.z};
    return mul_v;
}

vector add(vector v1, vector v2){
    vector added = {v1.x+v2.x, v1.y+v2.y,v1.z+v2.z};
    return added;
}

vector a = {0, 0, 0};
vector velocity = {0, 7560, 0};
vector position = {2e5,0,7e6};
int time;

void print_vector(vector v){
    printf("%.30f,%.30f,%.30f\n", v.x,v.y,v.z);
}

vector get_gravity(long double x, long double y, long double z){
    vector r = {x, y, z};
    long double r_mag = mag(r);
    
    long double ax = ((u*J2*(pow(R_e, 2)))/2)*((15*x*pow(z, 2))/pow(r_mag, 7) - (3*x/pow(r_mag, 5)));
    long double ay = ((u*J2*(pow(R_e, 2)))/2)*((15*y*pow(z, 2))/pow(r_mag, 7) - (3*y/pow(r_mag, 5)));
    long double az = ((u*J2*(pow(R_e, 2)))/2)*((15*pow(z, 3))/pow(r_mag, 7) - (9*z/pow(r_mag, 5)));

    vector a_j2 = {ax, ay, az};
    vector a = scalar_mult(-((G*M_e)/(pow(r_mag, 3))), r);
    return add(a_j2,a);
    //return a;
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
    return 2600000+time*dt;
}

long double get_time(int julian_date){
    return (julian_date - 2451545)/36525;
}

long double dot(vector v1, vector v2){
    return v1.x*v2.x + v1.y*v2.y + v1.z* v2.z;
}

vector get_SRP(vector position,long double time,vector* normals){
    long double phi = 280.460 + 36000.771*time;
    long double M = 357.5277233 + 35999.05034*time; 
    long double phi_ecl = phi + 1.914666471*sin(M) + 0.019994643*sin(2*M); 
    long double epsilon = 23.439291 - 0.0130042*time;

    vector earth_sun_unit = {cos(phi_ecl), cos(epsilon)*sin(phi_ecl), sin(epsilon)*sin(phi_ecl)};
    vector sat_sun = add(scalar_mult(dist_to_sun*epsilon,earth_sun_unit),scalar_mult(-1,position));

    vector sat_sun_unit = scalar_mult(1/mag(sat_sun),sat_sun);
    vector r_unit = scalar_mult(1/mag(position),position);

    long double theta = acos(dot(r_unit, earth_sun_unit));

    if ((theta > (pi/2)) && (mag(position)*sin(theta) < R_e)) {
        printf("FSRP dying\n");
        vector zero = {0, 0, 0};
        return zero;
    }

    long double P_weirdthing = S_weirdthing/(c*pow(mag(sat_sun), 2));

    long double costheta1 = fabs(dot(normals[0], sat_sun_unit));
    long double projArea1 = A1*costheta1;
    vector FSRP_face1 = scalar_mult(-P_weirdthing*projArea1, add(scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]), scalar_mult((1-R_spec), sat_sun_unit)));

//scalar_mult(-P_weirdthing*projArea1, add(scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]), scalar_mult((1-R_spec), sat_sun_unit)))
//add(scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]), scalar_mult((1-R_spec), sat_sun_unit))
//scalar_mult(2*(R_diff/3 + R_spec*costheta1), normals[0]);
//scalar_mult((1-R_spec), sat_sun_unit)


    long double costheta2 = fabs(dot(normals[1], sat_sun_unit));
    long double projArea2 = A2*costheta2;
    vector FSRP_face2 = scalar_mult(-P_weirdthing*projArea2, add(scalar_mult(2*(R_diff/3 + R_spec*costheta2), normals[1]), scalar_mult((1-R_spec), sat_sun_unit)));

    long double costheta3 = fabs(dot(normals[2], sat_sun_unit));
    long double projArea3 = A3*costheta3;
    vector FSRP_face3 = scalar_mult(-P_weirdthing*projArea3, add(scalar_mult(2*(R_diff/3 + R_spec*costheta3), normals[2]), scalar_mult((1-R_spec), sat_sun_unit)));

    vector FSRP = add(FSRP_face1, add(FSRP_face2, FSRP_face3));
    printf("FSRP FACES:\n");
    print_vector(scalar_mult(10e25, FSRP_face1));
    print_vector(scalar_mult(10e25, FSRP_face2));
    print_vector(scalar_mult(10e25, FSRP_face3));
    return FSRP;
}

vector drag(long double rho, vector* normals, vector velocity){
    vector unit_velocity = scalar_mult(1/mag(velocity), velocity);
    // print_vector(unit_velocity);

    long double costheta1 = fabs(dot(normals[0], unit_velocity));
    long double projArea1 = A1*costheta1;
    long double costheta2 = fabs(dot(normals[1], unit_velocity));
    long double projArea2 = A2*costheta2;
    long double costheta3 = fabs(dot(normals[2], unit_velocity));
    long double projArea3 = A3*costheta3;
    if(projArea3 != projArea3){
        exit(0);
    }
    long double A = projArea1 + projArea2 + projArea3;
    vector F_D = scalar_mult(((-rho*mag(velocity)*C_D*A)/2), velocity);
    // printf("rho %f\n", rho);
    return F_D;
}

long double get_rho(vector position){
    int h0[] = {
        0,25000,30000,35000,40000,45000,50000,55000,60000,65000,70000,75000,80000,85000,90000,95000,100000,110000,120000,130000,140000,150000,160000,180000,200000,250000,300000,350000,400000,450000,500000,600000,700000,800000,900000,1000000
    };

    long double rho0[] = {
        1.225,0.03899,0.01774,0.008279,0.003972,0.001995,0.001057,0.0005821,0.0003206,0.0001718,0.0000877,0.00004178,0.00001905,0.000008337,0.000003396,0.000001343,5.30E-07,9.66E-08,2.44E-08,8.48E-09,3.85E-09,2.07E-09,1.22E-09,5.46E-10,2.79E-10,7.25E-11,2.42E-11,9.16E-12,3.73E-12,1.59E-12,6.97E-13,1.45E-13,3.61E-14,1.17E-14,5.25E-15,3.02E-15
    };

    long double H[] = {
        8440,6490,6750,7070,7470,7830,7950,7730,7290,6810,6330,6000,5700,5410,5380,5740,6150,8060,11600,16100,20600,24600,26300,33200,38500,46900,52500,56400,59400,62200,65800,79000,109000,164000,225000,268000
    };

    long double h = mag(position) - R_e;
    int i = 0;
    for (i = 0; i<36; i++) {
        if (h<=h0[i]) break;
    }
    i--;

    if(i < 0){
        printf("Why life\n");
        return 0;
    }

    long double rho = rho0[i]*exp(-((h-(h0[i]))/H[i]));
    //printf("rho %Lf\n",rho0[i]);
    //printf("bruh %Lf\n",3.02e-15);
    //printf("h %Lf\n",h);
    //printf("exp %Lf",exp(-((h-(h0[i]))/H[i])));
    return rho;
}

int main(){
    FILE *fptr;
    fptr = fopen("Values.txt", "w");
    vector *pos = (vector*)malloc((total_time/dt)* sizeof(vector));
    for (time = 1; time<=total_time; time+=dt){
        //printf("mag %f\n", mag(position));
        vector a_t = get_gravity(position.x, position.y, position.z);
        vector *normals = (vector*)malloc(3* sizeof(vector));
        get_normals(normals);
        vector a_drag = drag(get_rho(position), normals, velocity);
        vector F_SRP = get_SRP(position, get_time(get_julian()), normals);
        vector a_SRP = scalar_mult(1/mass, F_SRP);
        a = add(a_t, add(a_drag, a_SRP));
        printf("ASRP: ");
        print_vector(scalar_mult(10e25,a_SRP));
        //a = add(a_t,a_drag);
        //a = a_t;
        velocity = add(velocity, scalar_mult(dt, a));
        position = add(position, scalar_mult(dt, velocity));
        pos[(time-1)/dt] = position; 
        
        fprintf(fptr, "%f %f %f\n", position.x, position.y, position.z);
    }
    fclose(fptr);
}