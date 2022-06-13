#include <iostream>
//#include <jakobi.h>

#include <Eigen/Dense>
#include <iostream>
#include <math.h>

#define pi (2 * acos(0.0))

using namespace std;

// DH matrika funkcija
Eigen::Matrix4d DH(double alpha, double a, double d, double theta)
{

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

Eigen::MatrixXd GeometricJakobian(double t0, double t1, double t2, double t3, double t4, double t5, double t6)
{

    Eigen::Matrix4d A1, A2, A3, A4, A5, A6, A7;
    Eigen::Matrix4d T2, T3, T4, T5, T6, T7;
    Eigen::Vector3d z0, z1, z2, z3, z4, z5, z6;
    Eigen::Vector3d p0, p1, p2, p3, p4, p5, p6;
    Eigen::Vector3d P;
    Eigen::MatrixXd jakobijeva_matrika(6, 7);
    Eigen::Vector3d cross1, cross2, cross3, cross4, cross5, cross6, cross7;

    A1 = DH(0, 0, 0.333, t0);
    A2 = DH(-pi / 2, 0, 0, t1);
    A3 = DH(pi / 2, 0, 0.316, t2);
    A4 = DH(pi / 2, 0.0825, 0, t3);
    A5 = DH(-pi / 2, -0.0825, 0.384, t4);
    A6 = DH(pi / 2, 0, 0, t5);
    A7 = DH(pi / 2, 0.088, 0.107, t6);

    T2 = A1 * A2;
    T3 = A1 * A2 * A3;
    T4 = A1 * A2 * A3 * A4;
    T5 = A1 * A2 * A3 * A4 * A5;
    T6 = A1 * A2 * A3 * A4 * A5 * A6;
    T7 = A1 * A2 * A3 * A4 * A5 * A6 * A7;

    z0 << 0, 0, 1;
    z1 << A1(0, 2), A1(1, 2), A1(2, 2);
    z2 << T2(0, 2), T2(1, 2), T2(2, 2);
    z3 << T3(0, 2), T3(1, 2), T3(2, 2);
    z4 << T4(0, 2), T4(1, 2), T4(2, 2);
    z5 << T5(0, 2), T5(1, 2), T5(2, 2);
    z6 << T6(0, 2), T6(1, 2), T6(2, 2);

    p0 << 0, 0, 0;
    p1 << A1(0, 3), A1(1, 3), A1(2, 3);
    p2 << T2(0, 3), T2(1, 3), T2(2, 3);
    p3 << T3(0, 3), T3(1, 3), T3(2, 3);
    p4 << T4(0, 3), T4(1, 3), T4(2, 3);
    p5 << T5(0, 3), T5(1, 3), T5(2, 3);
    p6 << T6(0, 3), T6(1, 3), T6(2, 3);

    P << T7(0, 3), T7(1, 3), T7(2, 3);

    cross1 << z0.cross(P - p0);
    cross2 << z1.cross(P - p1);
    cross3 << z2.cross(P - p2);
    cross4 << z3.cross(P - p3);
    cross5 << z4.cross(P - p4);
    cross6 << z5.cross(P - p5);
    cross7 << z6.cross(P - p6);

    jakobijeva_matrika << cross1(0), cross2(0), cross3(0), cross4(0), cross5(0), cross6(0), cross7(0),
                          cross1(1), cross2(1), cross3(1), cross4(1), cross5(1), cross6(1), cross7(1),
                          cross1(2), cross2(2), cross3(2), cross4(2), cross5(2), cross6(2), cross7(2),
                          z0(0),z1(0),z2(0),z3(0),z4(0),z5(0),z6(0),
                          z0(1),z1(1),z2(1),z3(1),z4(1),z5(1),z6(1),
                          z0(2),z1(2),z2(2),z3(2),z4(2),z5(2),z6(2);

    return jakobijeva_matrika;
}

int main(){
    std::cout<< GeometricJakobian(0,0,0,0,0,0,0) << std::endl;
}