#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#define pi (2 * acos(0.0))

using namespace std;

//DH matrika funkcija 
Eigen::Matrix4d DH(double a, double d, double alpha, double theta){

    Eigen::Matrix4d matrika;

    double ct = cos(theta);
    double st = sin(theta);
    double ca = cos(alpha);
    double sa = sin(alpha);

    matrika << ct, -st, 0, a,
               st * ca, ct * ca, -sa, -d * sa,
               st * sa, ct * sa, ca, d * ca,
                0, 0, 0, 1;

    return matrika;
}

Eigen::MatrixXd GeometricJakobian(double t0, double t1, double t2, double t3, double t4, double t5, double t6){

    Eigen::Matrix4d A1, A2, A3, A4, A5, A6, A7;
    Eigen::Matrix4d T2, T3, T4, T5, T6, T7;
    Eigen::Vector3d z0, z1 ,z2, z3, z4, z5, z6;
    Eigen::Vector3d p0, p1, p2, p3, p4, p5, p6;

    A1 = DH(0, 0, 0.333, t0);
    A2 = DH(-pi / 2, 0, 0, t1);
    A3 = DH(pi / 2, 0, 0.316, t2);
    A4 = DH(pi / 2, 0.0825, 0, t3);
    A5 = DH(-pi / 2, -0.0825, 0.384, t4);
    A6 = DH(pi / 2, 0, 0, t5);
    A7 = DH(pi / 2, 0.088, 0.107, t6);

    T2 = A1*A2;
    T3 = A1*A2*A3;
    T4 = A1*A2*A3*A4;
    T5 = A1*A2*A3*A4*A5;
    T6 = A1*A2*A3*A4*A5*A6;
    T7 = A1*A2*A3*A4*A5*A6*A7;

    z0 << 0,0,1;
    



}