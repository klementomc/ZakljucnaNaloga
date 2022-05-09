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

vector<double> seznam = {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0}; // nastavi prave vrednosti kotov

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
Eigen::MatrixXd DirektnaKinematika(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double theta7)
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

int main(){
    Eigen::MatrixXd MAT = DirektnaKinematika(0, 0, 0, 0, 0, 0, 0);
    cout << DirektnaKinematika(0, 0, 0, 0, 0, 0, 0) << endl;
    cout << DirektnaKinematika(0, 1.0323, 0, 0.8247, 0, 0.2076, 0) << endl;
    cout << RotMat2rpy(MAT) << endl;
    cout << DirektnaKinematika(0.512333, 0.632566, 0.653824, 0.797273, -0.687335, 0.106895, -0.359744) << endl;
}