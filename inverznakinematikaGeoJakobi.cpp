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

#define n_iter 1500
#define thershP 0.001
#define thresh0 0.001
#define koef_alpha 0.6
#define pi 2 * acos(0.0) // izračunamo koliko je Pi

using namespace std;

int stevec = 0;

typedef Eigen::Matrix<double, 3, 1> vektor3d; // typdef za seznam
typedef Eigen::Matrix<double, 7, 1> vektor7d; // typdef za seznam
// vector<double> seznam = {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0}; // nastavi prave vrednosti kotov

vector<double> seznam = {0, -17.2 * pi / 180, 0, -126 * pi / 180, 0, 115 * pi / 180, 45 * pi / 180}; // nastavi prave vrednosti kotov

// Izračun psevdoinverzne matrike Moor-Penrose inverse
// method for calculating the pseudo-Inverse as recommended by Eigen developers
template <typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
{
    Eigen::JacobiSVD<_Matrix_Type_> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

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
        z0(0), z1(0), z2(0), z3(0), z4(0), z5(0), z6(0),
        z0(1), z1(1), z2(1), z3(1), z4(1), z5(1), z6(1),
        z0(2), z1(2), z2(2), z3(2), z4(2), z5(2), z6(2);

    return jakobijeva_matrika;
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
    // cout << "1----------------------------" << endl;
    //  matrika KP in KO
    Eigen::MatrixXd Kp_Ko(6, 6);

    Kp_Ko << 3, 0, 0, 0, 0, 0,
        0, 3, 0, 0, 0, 0,
        0, 0, 3, 0, 0, 0,
        0, 0, 0, 3, 0, 0,
        0, 0, 0, 0, 3, 0,
        0, 0, 0, 0, 0, 3;

    // naredimo kvaternion vektor iz željenih kotov
    Eigen::Matrix3d Zacetna_rotacijska_matrika;
    Zacetna_rotacijska_matrika = RotMat(T_cilj(3), T_cilj(4), T_cilj(5));
    Eigen::Vector4d Zacetni_kvaternion_vektor = RotMat2Qvec(Zacetna_rotacijska_matrika);
    Eigen::Vector3d Zacetni_kvaternion_vektor_ijk;
    Zacetni_kvaternion_vektor_ijk << Zacetni_kvaternion_vektor(1), Zacetni_kvaternion_vektor(2), Zacetni_kvaternion_vektor(3);

    // cout<< "2----------------------------"<<endl;
    float neP = 1.00; // vrednosti za izračun napake
    float neO = 1.00; // vrednosti za izračun napake

    Eigen::Vector3d Zeljena_pozicija;
    Zeljena_pozicija << T_cilj(0), T_cilj(1), T_cilj(2);

    Eigen::Vector3d eP;
    Eigen::Vector3d eO;
    Eigen::Vector3d Trenutna_pozicija;
    Eigen::Vector3d Trenutna_orientacija;
    Eigen::Vector3d kvaternion_vektor_trenutni_ijk;
    Eigen::VectorXd dq(7);

    // cout << "3----------------------------" << endl;
    while (!(neP < thershP && neO < thresh0) && stevec < n_iter)
    {

        Eigen::MatrixXd J(6, 7);
        J = GeometricJakobian(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));
        Eigen::MatrixXd FK(4, 4);
        FK = DirektnaKinematika(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));
        // cout << "4----------------------------" << endl;

        // kvaternion vektor trenutnega stanja
        Eigen::Vector4d kvaternion_vektor_trenutni = RotMat2Qvec(FK);
        kvaternion_vektor_trenutni_ijk << kvaternion_vektor_trenutni(1), kvaternion_vektor_trenutni(2), kvaternion_vektor_trenutni(3);
        // cout << "5----------------------------" << endl;
        Trenutna_pozicija << FK(0, 3), FK(1, 3), FK(2, 3);

        eP = Zeljena_pozicija - Trenutna_pozicija;
        eO = -Zacetni_kvaternion_vektor(0) * kvaternion_vektor_trenutni_ijk + kvaternion_vektor_trenutni(0) * Zacetni_kvaternion_vektor_ijk - Zacetni_kvaternion_vektor_ijk.cross(kvaternion_vektor_trenutni_ijk);

        Eigen::VectorXd vektor_zdruzen_eP_eO(6);
        vektor_zdruzen_eP_eO << eP(0), eP(1), eP(2), eO(0), eO(1), eO(2);
        Eigen::MatrixXd Pseudo_J(7, 6); //
        Pseudo_J = pseudoInverse(J);
        dq = koef_alpha * (Pseudo_J * (Kp_Ko * vektor_zdruzen_eP_eO));

        TrenutniKoti = TrenutniKoti + dq;

        // cout << dq.transpose() << endl;

        // tu pride koda za prilagoitev sklepov

        for (int i = 0; i <= 6; i++)
        {
            if (TrenutniKoti(i) < -2 * pi)
            {
                while (TrenutniKoti(i) < 0)
                {
                    TrenutniKoti(i) = TrenutniKoti(i) + (2 * pi);
                }
                if (TrenutniKoti(i) > 0)
                {
                    TrenutniKoti(i) = TrenutniKoti(i) - (2 * pi);
                } // pogoj da se nebo predznak spremenu
            }
            if (TrenutniKoti(i) > 2 * pi)
            {
                while (TrenutniKoti(i) > 0)
                {
                    TrenutniKoti(i) = TrenutniKoti(i) - (2 * pi);
                }
                if (TrenutniKoti(i) < 0)
                {
                    TrenutniKoti(i) = TrenutniKoti(i) + (2 * pi);
                } // pogoj da se nebo predznak spremenu
            }
        }

        /*for (int i = 0; i <= 6; i++)
        {
            if (TrenutniKoti(i) < MIN[i])
            {
                TrenutniKoti(i) = MIN[i];
            }

            if (TrenutniKoti(i) > MAX[i])
            {
                TrenutniKoti(i) = MAX[i];
            }
        }*/

        neP = sqrt(pow(eP(0), 2) + pow(eP(1), 2) + pow(eP(2), 2));
        neO = sqrt(pow(eO(0), 2) + pow(eO(1), 2) + pow(eO(2), 2));

        stevec += 1;
    }
    cout << stevec << endl;
    stevec = 0;
    for (int u = 0; u < 7; u++)
    {
        seznam[u] = TrenutniKoti(u);
    }
    return TrenutniKoti;
}

int main()
{
    cout << "a" << endl;
    cout << "______________________________________________" << endl;
    cout << "" << endl;
    Eigen::VectorXd test(6);
    test << 0.6, 0.7, 0.3, pi, 0, pi;
    cout << InverznaKinematika(test) << endl;
    cout << "______________________________________________" << endl;
    cout << "" << endl;
    cout << DirektnaKinematika(seznam[0], seznam[1], seznam[2], seznam[3], seznam[4], seznam[5], seznam[6]) << endl;
}