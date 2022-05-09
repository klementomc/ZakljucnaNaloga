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

#define n_iter 1000
#define thershP 10 ^ -4
#define thresh0 10 ^ -4
#define koef_alpha 0.1

using namespace std;

int stevec = 0;

typedef Eigen::Matrix<double, 3, 1> vektor3d; // typdef za seznam
typedef Eigen::Matrix<double, 7, 1> vektor7d; // typdef za seznam

vector<double> seznam = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // nastavi prave vrednosti kotov

// Izračun psevdoinverzne matrike Moor-Penrose inverse
// method for calculating the pseudo-Inverse as recommended by Eigen developers
template <typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
{
    Eigen::JacobiSVD<_Matrix_Type_> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

double pi = 2 * acos(0.0); // izračunamo koliko je Pi

// ustvarimo objekt sklepa, ki ima funkcijo izračunat denavit hartenbergovo matriko transformacij
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
            st * ca, ct * ca, -sa, -d * sa,
            st * sa, ct * sa, ca, d * ca,
            0, 0, 0, 1;

        return DenavitHartenberg;
    }
};

// funkcija za izračun direktne kinematike
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

    J << 87999999999997 * (((-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * cos(theta4) - sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * sin(theta4) - sin(theta1) * sin(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * cos(theta4) / 643371375338703 + 53078138465443 * sin(theta1) * sin(theta2) * sin(theta4) / 643371375338703 - 383999999999957 * sin(theta1) * sin(theta2) * cos(theta4) / 999999999999888 - 315999999999994 * sin(theta1) * sin(theta2) / 999999999999981 - 53078138465443 * sin(theta1) * cos(theta2) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta3) * cos(theta1) / 643371375338703, 87999999999997 * ((-sin(theta2) * cos(theta1) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta1) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5) * cos(theta1)) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta1) * cos(theta3) + cos(theta1) * cos(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 + 383999999999957 * sin(theta2) * sin(theta4) * cos(theta1) * cos(theta3) / 999999999999888 + 53078138465443 * sin(theta2) * cos(theta1) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta2) * cos(theta1) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta4) * cos(theta1) * cos(theta2) / 643371375338703 + 383999999999957 * cos(theta1) * cos(theta2) * cos(theta4) / 999999999999888 + 315999999999994 * cos(theta1) * cos(theta2) / 999999999999981, 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta5) + (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta4) * cos(theta5)) * cos(theta6) / 999999999999966 - 87999999999997 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta3) * cos(theta1) * cos(theta2) / 643371375338703, 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) - sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta6) / 999999999999966 + 53078138465443 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) / 643371375338703 - 383999999999957 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta2) * sin(theta4) * cos(theta1) / 999999999999888 - 53078138465443 * sin(theta2) * cos(theta1) * cos(theta4) / 643371375338703, 87999999999997 * (-((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 87999999999997 * (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta3) / 643371375338703 - 53078138465443 * sin(theta2) * sin(theta4) * cos(theta1) / 643371375338703 + 383999999999957 * sin(theta2) * cos(theta1) * cos(theta4) / 999999999999888 + 315999999999994 * sin(theta2) * cos(theta1) / 999999999999981 + 53078138465443 * cos(theta1) * cos(theta2) * cos(theta3) / 643371375338703, 87999999999997 * ((-sin(theta1) * sin(theta2) * cos(theta3) * cos(theta4) + sin(theta1) * sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta1) * sin(theta2) * sin(theta3) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta1) * sin(theta2) * sin(theta4) * cos(theta3) + sin(theta1) * cos(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 + 383999999999957 * sin(theta1) * sin(theta2) * sin(theta4) * cos(theta3) / 999999999999888 + 53078138465443 * sin(theta1) * sin(theta2) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta2) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta4) * cos(theta2) / 643371375338703 + 383999999999957 * sin(theta1) * cos(theta2) * cos(theta4) / 999999999999888 + 315999999999994 * sin(theta1) * cos(theta2) / 999999999999981, 87999999999997 * ((-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * cos(theta4) * cos(theta5) - (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta5)) * cos(theta6) / 999999999999966 - 87999999999997 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta3) * cos(theta2) / 643371375338703 + 53078138465443 * cos(theta1) * cos(theta3) / 643371375338703, 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) - sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta6) / 999999999999966 + 53078138465443 * (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) / 643371375338703 - 383999999999957 * (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta1) * sin(theta2) * sin(theta4) / 999999999999888 - 53078138465443 * sin(theta1) * sin(theta2) * cos(theta4) / 643371375338703, 87999999999997 * (-((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * (((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 0, 87999999999997 * ((-sin(theta2) * sin(theta4) - cos(theta2) * cos(theta3) * cos(theta4)) * cos(theta5) + sin(theta3) * sin(theta5) * cos(theta2)) * cos(theta6) / 999999999999966 + 87999999999997 * (-sin(theta2) * cos(theta4) + sin(theta4) * cos(theta2) * cos(theta3)) * sin(theta6) / 999999999999966 + 53078138465443 * sin(theta2) * sin(theta4) / 643371375338703 - 383999999999957 * sin(theta2) * cos(theta4) / 999999999999888 - 315999999999994 * sin(theta2) / 999999999999981 + 383999999999957 * sin(theta4) * cos(theta2) * cos(theta3) / 999999999999888 + 53078138465443 * cos(theta2) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * cos(theta2) * cos(theta3) / 643371375338703, 87999999999997 * (sin(theta2) * sin(theta3) * cos(theta4) * cos(theta5) + sin(theta2) * sin(theta5) * cos(theta3)) * cos(theta6) / 999999999999966 - 87999999999997 * sin(theta2) * sin(theta3) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * sin(theta2) * sin(theta3) * sin(theta4) / 999999999999888 - 53078138465443 * sin(theta2) * sin(theta3) * cos(theta4) / 643371375338703 + 53078138465443 * sin(theta2) * sin(theta3) / 643371375338703, 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * cos(theta3) * cos(theta4) - sin(theta4) * cos(theta2)) * sin(theta6) / 999999999999966 - 53078138465443 * sin(theta2) * sin(theta4) * cos(theta3) / 643371375338703 + 383999999999957 * sin(theta2) * cos(theta3) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta4) * cos(theta2) / 999999999999888 - 53078138465443 * cos(theta2) * cos(theta4) / 643371375338703, 87999999999997 * (-(-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * sin(theta5) + sin(theta2) * sin(theta3) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * ((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 0, -sin(theta2) * cos(theta1), -sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2), -(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4), -((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta5), -(((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * sin(theta6) + (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6), -((((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * cos(theta6) + (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * sin(theta6)) * sin(theta7) + (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta5) + (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta5)) * cos(theta7), 0, -sin(theta1) * sin(theta2), -sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3), -(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4), -((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * cos(theta5), -(((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * sin(theta6) + (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6), -((((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * cos(theta6) + (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * sin(theta6)) * sin(theta7) + (((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta5) + (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * cos(theta5)) * cos(theta7), 1, -cos(theta2), sin(theta2) * sin(theta3), sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4), -(-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * sin(theta5) + sin(theta2) * sin(theta3) * cos(theta5), -((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * sin(theta6) + (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6), -(((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * cos(theta6) + (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * sin(theta6)) * sin(theta7) + ((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * sin(theta5) - sin(theta2) * sin(theta3) * cos(theta5)) * cos(theta7);

    return J;
}

// to ustvari kvaternion
Eigen::Quaterniond kvaternion(Eigen::MatrixXd M)
{

    Eigen::Matrix3d Nm;

    Nm << M(0, 0), M(0, 1), M(0, 2),
        M(1, 0), M(1, 1), M(1, 2),
        M(2, 0), M(2, 1), M(2, 2);

    Eigen::Quaterniond kvaternion(Nm);

    return kvaternion;
}

// naredi vektor iz rotacijse matrike oblike {n,eta}, kjer n predstavlja realni del eta pa imaginarni

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

Eigen::VectorXd InverznaKinematika(Eigen::VectorXd T_cilj)
{

    std::vector<double> MAX = {sklep1_MAX, sklep2_MAX, sklep3_MAX, sklep4_MAX, sklep5_MAX, sklep6_MAX, sklep7_MAX};
    std::vector<double> MIN = {sklep1_MIN, sklep2_MIN, sklep3_MIN, sklep4_MIN, sklep4_MIN, sklep6_MIN, sklep7_MIN};
    Eigen::VectorXd TrenutniKoti(7);
    for (int i = 0; i < 7; i++)
    {
        TrenutniKoti(i) = seznam[i];
    }
    //cout << "1----------------------------" << endl;
    // matrika KP in KO
    Eigen::MatrixXd Kp_Ko(6, 6);

    Kp_Ko << 1, 0, 0, 0, 0, 0,
             0, 1, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 0,
             0, 0, 0, 1, 0, 0,
             0, 0, 0, 0, 1, 0,
             0, 0, 0, 0, 0, 1;

    // naredimo kvaternion vektor iz željenih kotov
    Eigen::Matrix3d Zacetna_rotacijska_matrika;
    Zacetna_rotacijska_matrika = RotMat(T_cilj(3), T_cilj(4), T_cilj(5));
    Eigen::Vector4d Zacetni_kvaternion_vektor = RotMat2Qvec(Zacetna_rotacijska_matrika);
    Eigen::Vector3d Zacetni_kvaternion_vektor_ijk;
    Zacetni_kvaternion_vektor_ijk << Zacetni_kvaternion_vektor(1), Zacetni_kvaternion_vektor(2), Zacetni_kvaternion_vektor(3);

    //cout<< "2----------------------------"<<endl;
    float neP = 1.00; // vrednosti za izračun napake
    float neO = 1.00; // vrednosti za izračun napake

    Eigen::Vector3d Zeljena_pozicija;
    Zeljena_pozicija << T_cilj(0), T_cilj(1), T_cilj(2);

    Eigen::Vector3d eP;
    Eigen::Vector3d eO;
    Eigen::Vector3d Trenutna_pozicija;
    Eigen::Vector3d Trenutna_orientacija;
    Eigen::Vector3d kvaternion_vektor_trenutni_ijk;

    //cout << "3----------------------------" << endl;
    while (stevec < n_iter /*&& (neP < thershP && neO < thresh0)*/)
    {

        Eigen::MatrixXd J(6, 7);
        J = Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));
        Eigen::MatrixXd FK(4, 4);
        FK = DirektnaKinematika(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));
        //cout << "4----------------------------" << endl;

        // kvaternion vektor trenutnega stanja
        Eigen::Vector4d kvaternion_vektor_trenutni = RotMat2Qvec(FK);
        kvaternion_vektor_trenutni_ijk << kvaternion_vektor_trenutni(1), kvaternion_vektor_trenutni(2), kvaternion_vektor_trenutni(3);
        //cout << "5----------------------------" << endl;
        Trenutna_pozicija << FK(0, 3), FK(1, 3), FK(2, 3);

        eP = Zeljena_pozicija - Trenutna_pozicija;
        eO = -Zacetni_kvaternion_vektor(0) * kvaternion_vektor_trenutni_ijk + kvaternion_vektor_trenutni(0) * Zacetni_kvaternion_vektor_ijk - Zacetni_kvaternion_vektor_ijk.cross(kvaternion_vektor_trenutni_ijk);

        Eigen::VectorXd vektor_zdruzen_eP_eO(6);
        vektor_zdruzen_eP_eO << eP(0), eP(1), eP(2), eO(0), eO(1), eO(2);
        Eigen::MatrixXd Pseudo_J(7, 6); //
        Pseudo_J = pseudoInverse(J);
        TrenutniKoti = TrenutniKoti + koef_alpha * (Pseudo_J * (Kp_Ko * vektor_zdruzen_eP_eO));

        neP = eP.maxCoeff();
        neO = eO.maxCoeff();

        stevec += 1;
    }
    stevec = 0;
    for (int u = 0; u < 7; u++)
    {
        seznam[u] = TrenutniKoti(u);
    }
    return TrenutniKoti;
}

int main()
{
    cout << "a" <<endl;
    Eigen::VectorXd test(6);
    test << 0.4, 0, 0.71, 0, pi, 0;
    cout << InverznaKinematika(test);
}