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

using namespace std;

vector<double> seznam = { 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0 };//nastavi prave vrednosti kotov

template <typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_& a, double epsilon = std::numeric_limits<double>::epsilon())
{
    Eigen::JacobiSVD<_Matrix_Type_> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

double pi = 2 * acos(0.0); //izraèunamo koliko je Pi

// ustvarimo objekt sklepa, ki ima funkcijo izraèunat denavit hartenbergovo matriko transformacij
class Joint
{
public:
    double alpha;
    double a;
    double d;
    double theta;

    //konstruktor objekta

    Joint(double x, double y, double z, double xy) {
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

Eigen::MatrixXd DirektnaKinematika(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double theta7) {

    Eigen::MatrixXd DK(4, 4);// Matrika ki jo vrne funkcija

    //Pravilen vrstni red je alpha, a, d, theta
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

