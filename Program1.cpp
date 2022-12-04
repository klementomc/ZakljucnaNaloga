#include <iostream>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Core>
#include <math.h>
#include <vector>
#include <list>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <time.h>
#include <map>

#define sklep1_MAX 2.8973
#define sklep2_MAX 1.7628
#define sklep3_MAX 2.8973
#define sklep4_MAX -0.0698
#define sklep5_MAX 2.8973
#define sklep6_MAX 3.7525
#define sklep7_MAX 2.8973

#define sklep1_MIN -2.8973
#define sklep2_MIN -1.7628
#define sklep3_MIN -2.8973
#define sklep4_MIN -3.0718
#define sklep5_MIN -2.8973
#define sklep6_MIN -0.0175
#define sklep7_MIN -2.8973

#define rho 0.1
#define rho1 0.5
#define lambda 0.4

#define thershP 1e-3
#define threshO 1e-1
#define n_iter 1500
#define koef_alpha 0.03
#define nic 0

#define pi 2 * acos(0.0) 

using namespace std;

typedef Eigen::Matrix<double, 7, 1> vec7d;
typedef Eigen::Matrix<double, 3, 1> vec3d;
typedef Eigen::Matrix<double, 4, 1> vec4d;
typedef Eigen::Matrix<double, 6, 1> vec6d;

Eigen::MatrixXd pseudoInverse(Eigen::MatrixXd A)
{
    Eigen::MatrixXd AT = A.transpose();
    Eigen::MatrixXd M;
    Eigen::MatrixXd inverse = A * AT;
    Eigen::MatrixXd B; 
    B = inverse.inverse();

    M = AT * B;

    return M;
}

class Joint
{
public:
    double alpha;
    double a;
    double d;
    double theta;

    // konstruktor objekta

    Joint(double x, double y, double z, double xy)
    {
        alpha = x;
        a = y;
        d = z;
        theta = xy;
    }

    Eigen::MatrixXd DH()
    {

        Eigen::MatrixXd DenavitHartenberg(4, 4);

        double ct = cos(theta);
        double st = sin(theta);
        double ca = cos(alpha);
        double sa = sin(alpha);

        DenavitHartenberg << ct, -st, 0, a,
            st* ca, ct* ca, -sa, -d * sa,
            st* sa, ct* sa, ca, d* ca,
            0, 0, 0, 1;

        return DenavitHartenberg;
    }
};

Eigen::Matrix4d DirektnaKinematika(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double theta7)
{

    Eigen::MatrixXd DK(4, 4); // Matrika ki jo vrne funkcija

    // Pravilen vrstni red je alpha, a, d, theta
    Joint joint1 = Joint(0, 0, 0.333, theta1);
    Joint joint2 = Joint(-pi / 2, 0, 0, theta2);
    Joint joint3 = Joint(pi / 2, 0, 0.316, theta3);
    Joint joint4 = Joint(pi / 2, 0.0825, 0, theta4);
    Joint joint5 = Joint(-pi / 2, -0.0825, 0.384, theta5);
    Joint joint6 = Joint(pi / 2, 0, 0, theta6);
    Joint joint7 = Joint(pi / 2, 0.088, 0.107, theta7);

    DK = joint1.DH() * joint2.DH() * joint3.DH() * joint4.DH() * joint5.DH() * joint6.DH() * joint7.DH();

    return DK;
}

Eigen::MatrixXd Jakobi(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double theta7)
{
    Eigen::MatrixXd J(6, 7);

    J << 87999999999997 * (((-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * cos(theta4) - sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * sin(theta4) - sin(theta1) * sin(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * cos(theta4) / 643371375338703 + 53078138465443 * sin(theta1) * sin(theta2) * sin(theta4) / 643371375338703 - 383999999999957 * sin(theta1) * sin(theta2) * cos(theta4) / 999999999999888 - 315999999999994 * sin(theta1) * sin(theta2) / 999999999999981 - 53078138465443 * sin(theta1) * cos(theta2) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta3) * cos(theta1) / 643371375338703, 87999999999997 * ((-sin(theta2) * cos(theta1) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta1) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5) * cos(theta1)) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta1) * cos(theta3) + cos(theta1) * cos(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 + 383999999999957 * sin(theta2) * sin(theta4) * cos(theta1) * cos(theta3) / 999999999999888 + 53078138465443 * sin(theta2) * cos(theta1) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta2) * cos(theta1) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta4) * cos(theta1) * cos(theta2) / 643371375338703 + 383999999999957 * cos(theta1) * cos(theta2) * cos(theta4) / 999999999999888 + 315999999999994 * cos(theta1) * cos(theta2) / 999999999999981, 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta5) + (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta4) * cos(theta5)) * cos(theta6) / 999999999999966 - 87999999999997 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta3) * cos(theta1) * cos(theta2) / 643371375338703, 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) - sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta6) / 999999999999966 + 53078138465443 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) / 643371375338703 - 383999999999957 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta2) * sin(theta4) * cos(theta1) / 999999999999888 - 53078138465443 * sin(theta2) * cos(theta1) * cos(theta4) / 643371375338703, 87999999999997 * (-((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 87999999999997 * (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta3) / 643371375338703 - 53078138465443 * sin(theta2) * sin(theta4) * cos(theta1) / 643371375338703 + 383999999999957 * sin(theta2) * cos(theta1) * cos(theta4) / 999999999999888 + 315999999999994 * sin(theta2) * cos(theta1) / 999999999999981 + 53078138465443 * cos(theta1) * cos(theta2) * cos(theta3) / 643371375338703, 87999999999997 * ((-sin(theta1) * sin(theta2) * cos(theta3) * cos(theta4) + sin(theta1) * sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta1) * sin(theta2) * sin(theta3) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta1) * sin(theta2) * sin(theta4) * cos(theta3) + sin(theta1) * cos(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 + 383999999999957 * sin(theta1) * sin(theta2) * sin(theta4) * cos(theta3) / 999999999999888 + 53078138465443 * sin(theta1) * sin(theta2) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta2) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta4) * cos(theta2) / 643371375338703 + 383999999999957 * sin(theta1) * cos(theta2) * cos(theta4) / 999999999999888 + 315999999999994 * sin(theta1) * cos(theta2) / 999999999999981, 87999999999997 * ((-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * cos(theta4) * cos(theta5) - (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta5)) * cos(theta6) / 999999999999966 - 87999999999997 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta3) * cos(theta2) / 643371375338703 + 53078138465443 * cos(theta1) * cos(theta3) / 643371375338703, 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) - sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta6) / 999999999999966 + 53078138465443 * (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) / 643371375338703 - 383999999999957 * (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta1) * sin(theta2) * sin(theta4) / 999999999999888 - 53078138465443 * sin(theta1) * sin(theta2) * cos(theta4) / 643371375338703, 87999999999997 * (-((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * (((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 0, 87999999999997 * ((-sin(theta2) * sin(theta4) - cos(theta2) * cos(theta3) * cos(theta4)) * cos(theta5) + sin(theta3) * sin(theta5) * cos(theta2)) * cos(theta6) / 999999999999966 + 87999999999997 * (-sin(theta2) * cos(theta4) + sin(theta4) * cos(theta2) * cos(theta3)) * sin(theta6) / 999999999999966 + 53078138465443 * sin(theta2) * sin(theta4) / 643371375338703 - 383999999999957 * sin(theta2) * cos(theta4) / 999999999999888 - 315999999999994 * sin(theta2) / 999999999999981 + 383999999999957 * sin(theta4) * cos(theta2) * cos(theta3) / 999999999999888 + 53078138465443 * cos(theta2) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * cos(theta2) * cos(theta3) / 643371375338703, 87999999999997 * (sin(theta2) * sin(theta3) * cos(theta4) * cos(theta5) + sin(theta2) * sin(theta5) * cos(theta3)) * cos(theta6) / 999999999999966 - 87999999999997 * sin(theta2) * sin(theta3) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * sin(theta2) * sin(theta3) * sin(theta4) / 999999999999888 - 53078138465443 * sin(theta2) * sin(theta3) * cos(theta4) / 643371375338703 + 53078138465443 * sin(theta2) * sin(theta3) / 643371375338703, 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * cos(theta3) * cos(theta4) - sin(theta4) * cos(theta2)) * sin(theta6) / 999999999999966 - 53078138465443 * sin(theta2) * sin(theta4) * cos(theta3) / 643371375338703 + 383999999999957 * sin(theta2) * cos(theta3) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta4) * cos(theta2) / 999999999999888 - 53078138465443 * cos(theta2) * cos(theta4) / 643371375338703, 87999999999997 * (-(-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * sin(theta5) + sin(theta2) * sin(theta3) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * ((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 0, -sin(theta1), sin(theta2)* cos(theta1), sin(theta1)* cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2), -(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4), ((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1))* sin(theta5) + (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta5), (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5))* sin(theta6) - (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6), 0, cos(theta1), sin(theta1)* sin(theta2), sin(theta1)* sin(theta3)* cos(theta2) - cos(theta1) * cos(theta3), -(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4), ((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4))* sin(theta5) + (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * cos(theta5), (((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5))* sin(theta6) - (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6), 1, 0, cos(theta2), -sin(theta2) * sin(theta3), sin(theta2)* sin(theta4)* cos(theta3) + cos(theta2) * cos(theta4), (-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2))* sin(theta5) - sin(theta2) * sin(theta3) * cos(theta5), ((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5))* sin(theta6) - (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6);

    return J;
}

Eigen::Vector4d RotMat2Qvec(Eigen::MatrixXd M)
{
    Eigen::Matrix3d RotMat;
    RotMat << M(0, 0), M(0, 1), M(0, 2),
        M(1, 0), M(1, 1), M(1, 2),
        M(2, 0), M(2, 1), M(2, 2);

    Eigen::Quaterniond q(RotMat);
    Eigen::Vector4d vektor;
    vektor << q.w(), q.x(), q.y(), q.z();

    return vektor;
}

Eigen::Vector3d RotMat2rpy(Eigen::MatrixXd M)
{

    Eigen::Matrix3d RotMat;
    RotMat << M(0, 0), M(0, 1), M(0, 2),
        M(1, 0), M(1, 1), M(1, 2),
        M(2, 0), M(2, 1), M(2, 2);

    Eigen::Quaterniond q(RotMat);

    double roll = atan2(2 * q.y() * q.w() - 2 * q.x() * q.z(), 1 - 2 * q.y() * q.y() - 2 * q.z() * q.z());
    double pitch = atan2(2 * q.x() * q.w() - 2 * q.y() * q.z(), 1 - 2 * q.x() * q.x() - 2 * q.z() * q.z());
    double yaw = asin(2 * q.x() * q.y() + 2 * q.z() * q.w());

    Eigen::Vector3d rpy;
    rpy << roll, pitch, yaw;

    return rpy;
}

Eigen::Matrix3d RotMat(double x, double y, double z)
{
    Eigen::Matrix3d Rx;
    Eigen::Matrix3d Ry;
    Eigen::Matrix3d Rz;
    Eigen::Matrix3d R;

    Rx << 1, 0, 0,
        0, cos(x), -sin(x),
        0, sin(x), cos(x);

    Ry << cos(y), 0, sin(y),
        0, 1, 0 - sin(y), 0, cos(y);

    Rz << cos(z), -sin(z), 0,
        sin(z), cos(z), 0,
        0, 0, 1;

    R = Rz * Ry * Rx;

    return R;
}

double  kot_med_rotacijskima_matrikama(Eigen::Matrix3d A, Eigen::Matrix3d B) {

    Eigen::Matrix3d T;
    T = B.transpose();

    Eigen::Matrix3d R = A * T;

    //double trace = R(0, 0) + R(1, 1) + R(3, 3);

    double theta = 0.5;//sin((trace - 1) / 2);

    return theta;
}

Eigen::Matrix3d matrika_3x3(Eigen::MatrixXd M) {
    Eigen::Matrix3d RotMat;
    RotMat << M(0, 0), M(0, 1), M(0, 2),
        M(1, 0), M(1, 1), M(1, 2),
        M(2, 0), M(2, 1), M(2, 2);
    return RotMat;
}

// podaj kot w, x,y,z
Eigen::Matrix3d quat2RotMat(vec4d Q) {
    double q0 = Q(0);
    double q1 = Q(1);
    double q2 = Q(2);
    double q3 = Q(3);


    double r00 = 2 * (q0 * q0 + q1 * q1) - 1;
    double r01 = 2 * (q1 * q2 - q0 * q3);
    double r02 = 2 * (q1 * q3 + q0 * q2);


    double r10 = 2 * (q1 * q2 + q0 * q3);
    double r11 = 2 * (q0 * q0 + q2 * q2) - 1;
    double r12 = 2 * (q2 * q3 - q0 * q1);


    double r20 = 2 * (q1 * q3 - q0 * q2);
    double r21 = 2 * (q2 * q3 + q0 * q1);
    double r22 = 2 * (q0 * q0 + q3 * q3) - 1;

    Eigen::Matrix3d RotMat;
    RotMat << r00, r01, r02,
        r10, r11, r12,
        r20, r21, r22;

    return RotMat;
}

Eigen::Matrix4d DH(double alpha, double a, double d, double theta)
{

    Eigen::Matrix4d matrika;

    double ct = cos(theta);
    double st = sin(theta);
    double ca = cos(alpha);
    double sa = sin(alpha);

    matrika << ct, -st, 0, a,
        st* ca, ct* ca, -sa, -d * sa,
        st* sa, ct* sa, ca, d* ca,
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
        z0(0), z1(0), z2(0), z3(0), z4(0), z5(0), z6(0),
        z0(1), z1(1), z2(1), z3(1), z4(1), z5(1), z6(1),
        z0(2), z1(2), z2(2), z3(2), z4(2), z5(2), z6(2);

    return jakobijeva_matrika;
}

std::vector<double> MAX = { sklep1_MAX, sklep2_MAX, sklep3_MAX, sklep4_MAX, sklep5_MAX, sklep6_MAX, sklep7_MAX };
std::vector<double> MIN = { sklep1_MIN, sklep2_MIN, sklep3_MIN, sklep4_MIN, sklep4_MIN, sklep6_MIN, sklep7_MIN };

vec7d activation_and_sign_function(vec7d q_i, std::vector<double> q_min, std::vector<double> q_max) {

    std::vector<int> g_i;

    for (int i = 0; i <= 6; i++) {
        if (q_i(i) < q_min[i]) { g_i.push_back(-1); }
        if (q_max[i] < q_i(i)) { g_i.push_back(1); }
        else
        {
            g_i.push_back(0);
        }
    }

    vec7d g;
    g << g_i[0], g_i[1], g_i[2], g_i[3], g_i[4], g_i[5], g_i[6];

    return g;
}

double sigmoid_function_max(double q_i, double q_l0_max, double q_l1_min) {

    double lambda_max_i = 1 / (1 + exp(-12 * (q_i - q_l0_max) / (q_l1_min - q_l0_max) + 6));

    return lambda_max_i;
}

double sigmoid_function_min(double q_i, double q_l0_min, double q_l1_min) {

    double lambda_min_i = 1 / (1 + exp(-12 * (q_i - q_l0_min) / (q_l1_min - q_l0_min) + 6));

    return lambda_min_i;
}

Eigen::MatrixXd Pe(Eigen::MatrixXd J) {
    Eigen::MatrixXd In(7,7);
    In << 1, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 1; 0,
        0, 0, 0, 0, 0, 0, 1;

    Eigen::MatrixXd JT(7,6);
    JT = pseudoInverse(J);

    Eigen::MatrixXd A(7,7);
    A = JT * J;
    //cout << "---------------------------" << endl;
    //cout << A << endl;

    Eigen::MatrixXd B(7,7);
    B = In - A;
    return B;
}

Eigen::MatrixXd inverznaKinematika(vec3d pozicija, vec4d Zeljena_orientacija, vec7d ZacetniKoti) {

    int stevec = 0;

    vec7d TrenutniKoti = ZacetniKoti;
    vec3d eP;
    vec3d eO;
    vec7d q1;
    vec3d Trenutna_pozicija;
    vec3d Trenutna_orientacija;
    vec4d kvaternion_vektor_trenutni;
    vec3d kvaternion_vektor_trenutni_ijk;
    vec3d kvaternion_vektor_zacetni_ijk;

    kvaternion_vektor_zacetni_ijk << Zeljena_orientacija(1), Zeljena_orientacija(2), Zeljena_orientacija(3);

    double neP;
    double neO;
 

    std::vector<double> q_l0_min;
    std::vector<double> q_l0_max;
    std::vector<double> q_l1_min;
    std::vector<double> q_l1_max;
    std::vector<double> delta_q;

    for (int i = 0; i <= 6; i++) {
        double delta_q_i = MAX[i] - MIN[i];
        delta_q.push_back(delta_q_i);
    }

    for (int i = 0; i <= 6; i++) {

        double q_l0_min_i = MIN[i] + rho * delta_q[i];
        double q_l0_max_i = MAX[i] - rho * delta_q[i];

        q_l0_min.push_back(q_l0_min_i);
        q_l0_max.push_back(q_l0_max_i);
    }

    for (int i = 0; i <= 6; i++) {

        double q_l1_min_i = q_l0_min[i] - rho * rho1 * delta_q[i];
        double q_l1_max_i = q_l0_max[i] + rho * rho1 * delta_q[i];

        q_l1_min.push_back(q_l1_min_i);
        q_l1_max.push_back(q_l1_max_i);
    }


    Eigen::MatrixXd Kp_Ko(6, 6);
    Kp_Ko << 3, 0, 0, 0, 0, 0,
        0, 3, 0, 0, 0, 0,
        0, 0, 3, 0, 0, 0,
        0, 0, 0, 3, 0, 0,
        0, 0, 0, 0, 3, 0,
        0, 0, 0, 0, 0, 3;

    std::vector<double> napakaP;
    std::vector<double> napakaO;

    while (stevec < n_iter)
    {
        std::vector<double> lambda_l;

        //resetirej vektor ko prides do konca while loopa
        for (int i = 0; i <= 6; i++) {

            if (TrenutniKoti(i) < q_l1_min[i] || q_l1_max[i] < TrenutniKoti(i)) {
                lambda_l.push_back(i);
            }

            if (q_l1_min[i] <= TrenutniKoti(i) <= q_l0_min[i]) {
                double lambda_l_i = sigmoid_function_min(TrenutniKoti(i), q_l0_min[i], q_l1_min[i]);

                lambda_l.push_back(lambda_l_i);
            }

            if (q_l0_max[i] <= TrenutniKoti(i) <= q_l1_max[i]) {

                double lambda_l_i = sigmoid_function_max(TrenutniKoti(i), q_l0_max[i], q_l1_min[i]);

                lambda_l.push_back(lambda_l_i);
            }

            else
            {
                lambda_l.push_back(0);
            }
        }



        Eigen::MatrixXd J(6, 7);
        J = Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));
        Eigen::MatrixXd FK(4, 4);
        FK = DirektnaKinematika(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));

        Trenutna_pozicija << FK(0, 3), FK(1, 3), FK(2, 3);

        kvaternion_vektor_trenutni = RotMat2Qvec(FK);
        kvaternion_vektor_trenutni_ijk << kvaternion_vektor_trenutni(1), kvaternion_vektor_trenutni(2), kvaternion_vektor_trenutni(3);

        eP = pozicija - Trenutna_pozicija;
        eO = -Zeljena_orientacija(0) * kvaternion_vektor_trenutni_ijk + kvaternion_vektor_trenutni(0) * kvaternion_vektor_zacetni_ijk - kvaternion_vektor_zacetni_ijk.cross(kvaternion_vektor_trenutni_ijk);
        neP = sqrt(pow(eP(0), 2) + pow(eP(1), 2) + pow(eP(2), 2));
        std::cout << neP << endl;
        Eigen::Matrix3d Zeljena = quat2RotMat(Zeljena_orientacija);
        Eigen::Matrix3d Trenutna = matrika_3x3(FK);
        neO = sqrt(pow(eO(0), 2) + pow(eO(1), 2) + pow(eO(2), 2));
        napakaP.push_back(neP);
        napakaO.push_back(neO);

        //if (abs(neP) < thershP && abs(neO) < threshO){break;}

        //if (abs(neP) < thershP) { break; }

        vec6d vektor_zdruzen_eP_eO(6);
        vektor_zdruzen_eP_eO << eP, eO;

        //cout << vektor_zdruzen_eP_eO << endl;

        Eigen::MatrixXd Pseudo_J;
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(Kp_Ko * J);
         
        q1 = koef_alpha * svd.solve(Kp_Ko *vektor_zdruzen_eP_eO);
        //std::cout << q1 << endl;
        //std::cout << "--------" << endl;
        std::vector<double> q2;
        vec7d predznak;
        predznak = activation_and_sign_function(TrenutniKoti, q_l0_min, q_l0_max);
        //std::cout << predznak << endl;
        //std::cout << "-------------" << endl;
        vec6d predznak_brez_q7;
        predznak_brez_q7 << predznak(0), predznak(1), predznak(2), predznak(3), predznak(4), predznak(5);

        //cout << Pe(J) << endl;

        for (int i = 0; i <= 6; i++) {
            if (TrenutniKoti(i) < q_l1_min[i] || q_l1_max[i] < TrenutniKoti(i)) {
                double q2_i = (-(1 + lambda) * (abs(q1(i))) / abs((Pe(J) * predznak)(i)) * Pe(J) * predznak)(i);
                q2.push_back(q2_i);
            }
            if ((q_l1_min[i] <= TrenutniKoti(i) && TrenutniKoti(i) <= q_l0_min[i]) || (q_l0_max[i] <= TrenutniKoti(i) && TrenutniKoti(i) <= q_l1_max[i])){
                double q1_i = abs(q1(i));
                //std::cout << q1_i << endl;
                //std::cout << "--------" << endl;
                vec7d P_G = Pe(J) * predznak;
                //std::cout << P_G << endl;
                //std::cout << "--------" << endl;
                double P_G_i = abs(P_G(i));
                //std::cout << P_G_i << endl;
                //std::cout << "--------" << endl;
                double ulomek = q1_i / P_G_i;
                //std::cout << ulomek << endl;
                //std::cout << "--------" << endl;
                vec7d q2_vektor = -(1 + lambda) * ulomek * P_G;
                //std::cout << q2_vektor(i) << endl;
                q2.push_back(q2_vektor(i));
            }
            if(q_l0_min[i]<TrenutniKoti(i) || TrenutniKoti(i)<q_l0_max[i])
            {
                q2.push_back(nic);
            }
        }

        //cout << q2.size() << endl;

        vec7d eigen_q2;
        for (int i = 0; i <= 6; i++) {
            eigen_q2(i) = q2[i];
        }

        //std::cout << eigen_q2 << endl;

        vec7d dq;
        dq = q1 + eigen_q2;
        //cout << dq << endl;
        TrenutniKoti = TrenutniKoti + dq;
        //cout << TrenutniKoti << endl;

        for (int i = 0; i <= 6; i++) {
            if (TrenutniKoti(i) < -2 * pi) {
                while (TrenutniKoti(i) < 0) { TrenutniKoti(i) = TrenutniKoti(i) + (2 * pi); }
                if (TrenutniKoti(i) > 0) { TrenutniKoti(i) = TrenutniKoti(i) - (2 * pi); }//pogoj da se nebo predznak spremenu
            }
            if (TrenutniKoti(i) > 2 * pi) {
                while (TrenutniKoti(i) > 0) { TrenutniKoti(i) = TrenutniKoti(i) - (2 * pi); }
                if (TrenutniKoti(i) < 0) { TrenutniKoti(i) = TrenutniKoti(i) + (2 * pi); } // pogoj da se nebo predznak spremenu
            }
        }
        lambda_l.clear();
        q2.clear();

        stevec += 1;

    }


    Eigen::MatrixXd J_Transpose(7, 6);
    J_Transpose = Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6)).transpose();
    Eigen::MatrixXd A(7, 7);
    A = Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6)) * J_Transpose;

    double Yoshikawa = A.determinant();
    cout << "gibljivost -----------" << endl;
    cout << Yoshikawa << endl;
    cout << "----------------------" << endl;

    //TrenutniKoti(6) = 0.0497;
    std::cout << stevec << endl;
    std::cout << "_______________________________" << endl;
    std::cout << DirektnaKinematika(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6)) << endl;
    std::cout << "kvaternion-----------------------------" << endl;
    std::cout << RotMat2Qvec(DirektnaKinematika(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6))) << endl;
    std::cout << "koti sklepov-------------------------------" << endl;
    //std::cout << TrenutniKoti(0)<<", " << TrenutniKoti(1) << ", " << TrenutniKoti(2) << ", " << TrenutniKoti(3) << ", " << TrenutniKoti(4) << ", " << TrenutniKoti(5) << ", " << TrenutniKoti(6) << endl;

    return TrenutniKoti;
}



int main() {
    vec3d pozicija;
    vec4d Zeljena_orientacija;
    vec7d ZacetniKoti;

    pozicija << -0.063, -0.360, 0.346;
    Zeljena_orientacija << 0.607, -0.286, 0.741, -0.036;
    ZacetniKoti << -1.62, 0.039, -0.33, -1.85, 1.6, 1.05, -1.04;


    //std::cout << "a" << endl;
    std::cout << "______________________________________________" << endl;
    std::cout << "" << endl;
    std::cout << inverznaKinematika(pozicija, Zeljena_orientacija, ZacetniKoti) << endl;

}





