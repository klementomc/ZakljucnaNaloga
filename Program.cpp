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

vector<double> seznam = {0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0};//nastavi prave vrednosti kotov

// Izračun psevdoinverzne matrike Moor-Penrose inverse
// method for calculating the pseudo-Inverse as recommended by Eigen developers
template <typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
{
    Eigen::JacobiSVD<_Matrix_Type_> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

double pi = 2 * acos(0.0); //izračunamo koliko je Pi

// ustvarimo objekt sklepa, ki ima funkcijo izračunat denavit hartenbergovo matriko transformacij
class Joint
{
public:
    double alpha;
    double a;
    double d;
    double theta;

    //konstruktor objekta

    Joint(double x, double y, double z, double xy){
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

    DenavitHartenberg << ct, -st , 0, a,
                         st*ca, ct * ca, -sa, -d * sa,
                         st*sa, ct*sa, ca, d*ca,
                         0, 0, 0, 1;

    return DenavitHartenberg;
    }
};

//funkcija za izračun direktne kinematike
Eigen::MatrixXd DirektnaKinematika(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double theta7){

    Eigen::MatrixXd DK(4,4);// Matrika ki jo vrne funkcija

    //Pravilen vrstni red je alpha, a, d, theta
    Joint joint1 = Joint(0,0,0.333,theta1);
    Joint joint2 = Joint(-pi/2,0,0,theta2);
    Joint joint3 = Joint(pi/2,0,0.316,theta3);
    Joint joint4 = Joint(pi/2,0.0825,0,theta4);
    Joint joint5 = Joint(-pi/2,-0.0825,0.384,theta5);
    Joint joint6 = Joint(pi/2,0,0,theta6);
    Joint joint7 = Joint(pi/2,0.088,0.107,theta7);

    DK = joint1.DH() * joint2.DH() * joint3.DH() * joint4.DH() * joint5.DH() * joint6.DH() * joint7.DH();

    return DK;
}

Eigen::MatrixXd Jakobi(double theta1, double theta2, double theta3, double theta4, double theta5, double theta6, double theta7){

    Eigen::MatrixXd J(6,7);

    J << 87999999999997 * (((-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * cos(theta4) - sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * sin(theta4) - sin(theta1) * sin(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * cos(theta2) * cos(theta3) - sin(theta3) * cos(theta1)) * cos(theta4) / 643371375338703 + 53078138465443 * sin(theta1) * sin(theta2) * sin(theta4) / 643371375338703 - 383999999999957 * sin(theta1) * sin(theta2) * cos(theta4) / 999999999999888 - 315999999999994 * sin(theta1) * sin(theta2) / 999999999999981 - 53078138465443 * sin(theta1) * cos(theta2) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta3) * cos(theta1) / 643371375338703, 87999999999997 * ((-sin(theta2) * cos(theta1) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta1) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5) * cos(theta1)) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta1) * cos(theta3) + cos(theta1) * cos(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 + 383999999999957 * sin(theta2) * sin(theta4) * cos(theta1) * cos(theta3) / 999999999999888 + 53078138465443 * sin(theta2) * cos(theta1) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta2) * cos(theta1) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta4) * cos(theta1) * cos(theta2) / 643371375338703 + 383999999999957 * cos(theta1) * cos(theta2) * cos(theta4) / 999999999999888 + 315999999999994 * cos(theta1) * cos(theta2) / 999999999999981, 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta5) + (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta4) * cos(theta5)) * cos(theta6) / 999999999999966 - 87999999999997 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * cos(theta3) - sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta3) * cos(theta1) * cos(theta2) / 643371375338703, 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) - sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta6) / 999999999999966 + 53078138465443 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) / 643371375338703 - 383999999999957 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta2) * sin(theta4) * cos(theta1) / 999999999999888 - 53078138465443 * sin(theta2) * cos(theta1) * cos(theta4) / 643371375338703, 87999999999997 * (-((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * sin(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 87999999999997 * (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta3) / 643371375338703 - 53078138465443 * sin(theta2) * sin(theta4) * cos(theta1) / 643371375338703 + 383999999999957 * sin(theta2) * cos(theta1) * cos(theta4) / 999999999999888 + 315999999999994 * sin(theta2) * cos(theta1) / 999999999999981 + 53078138465443 * cos(theta1) * cos(theta2) * cos(theta3) / 643371375338703, 87999999999997 * ((-sin(theta1) * sin(theta2) * cos(theta3) * cos(theta4) + sin(theta1) * sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta1) * sin(theta2) * sin(theta3) * sin(theta5)) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta1) * sin(theta2) * sin(theta4) * cos(theta3) + sin(theta1) * cos(theta2) * cos(theta4)) * sin(theta6) / 999999999999966 + 383999999999957 * sin(theta1) * sin(theta2) * sin(theta4) * cos(theta3) / 999999999999888 + 53078138465443 * sin(theta1) * sin(theta2) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta2) * cos(theta3) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta4) * cos(theta2) / 643371375338703 + 383999999999957 * sin(theta1) * cos(theta2) * cos(theta4) / 999999999999888 + 315999999999994 * sin(theta1) * cos(theta2) / 999999999999981, 87999999999997 * ((-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * cos(theta4) * cos(theta5) - (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta5)) * cos(theta6) / 999999999999966 - 87999999999997 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * sin(theta4) / 999999999999888 - 53078138465443 * (-sin(theta1) * sin(theta3) * cos(theta2) + cos(theta1) * cos(theta3)) * cos(theta4) / 643371375338703 - 53078138465443 * sin(theta1) * sin(theta3) * cos(theta2) / 643371375338703 + 53078138465443 * cos(theta1) * cos(theta3) / 643371375338703, 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) - sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta6) / 999999999999966 + 53078138465443 * (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) / 643371375338703 - 383999999999957 * (sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta1) * sin(theta2) * sin(theta4) / 999999999999888 - 53078138465443 * sin(theta1) * sin(theta2) * cos(theta4) / 643371375338703, 87999999999997 * (-((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * sin(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * (((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 0, 87999999999997 * ((-sin(theta2) * sin(theta4) - cos(theta2) * cos(theta3) * cos(theta4)) * cos(theta5) + sin(theta3) * sin(theta5) * cos(theta2)) * cos(theta6) / 999999999999966 + 87999999999997 * (-sin(theta2) * cos(theta4) + sin(theta4) * cos(theta2) * cos(theta3)) * sin(theta6) / 999999999999966 + 53078138465443 * sin(theta2) * sin(theta4) / 643371375338703 - 383999999999957 * sin(theta2) * cos(theta4) / 999999999999888 - 315999999999994 * sin(theta2) / 999999999999981 + 383999999999957 * sin(theta4) * cos(theta2) * cos(theta3) / 999999999999888 + 53078138465443 * cos(theta2) * cos(theta3) * cos(theta4) / 643371375338703 - 53078138465443 * cos(theta2) * cos(theta3) / 643371375338703, 87999999999997 * (sin(theta2) * sin(theta3) * cos(theta4) * cos(theta5) + sin(theta2) * sin(theta5) * cos(theta3)) * cos(theta6) / 999999999999966 - 87999999999997 * sin(theta2) * sin(theta3) * sin(theta4) * sin(theta6) / 999999999999966 - 383999999999957 * sin(theta2) * sin(theta3) * sin(theta4) / 999999999999888 - 53078138465443 * sin(theta2) * sin(theta3) * cos(theta4) / 643371375338703 + 53078138465443 * sin(theta2) * sin(theta3) / 643371375338703, 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta5) * cos(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * cos(theta3) * cos(theta4) - sin(theta4) * cos(theta2)) * sin(theta6) / 999999999999966 - 53078138465443 * sin(theta2) * sin(theta4) * cos(theta3) / 643371375338703 + 383999999999957 * sin(theta2) * cos(theta3) * cos(theta4) / 999999999999888 - 383999999999957 * sin(theta4) * cos(theta2) / 999999999999888 - 53078138465443 * cos(theta2) * cos(theta4) / 643371375338703, 87999999999997 * (-(-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * sin(theta5) + sin(theta2) * sin(theta3) * cos(theta5)) * cos(theta6) / 999999999999966, -87999999999997 * ((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * sin(theta6) / 999999999999966 + 87999999999997 * (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6) / 999999999999966, 0, 0, -sin(theta2) * cos(theta1), sin(theta2) * cos(theta1), -(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4), -(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4), -(((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * sin(theta6) + (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6), (((-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * cos(theta4) + sin(theta2) * sin(theta4) * cos(theta1)) * cos(theta5) - (sin(theta1) * cos(theta3) + sin(theta3) * cos(theta1) * cos(theta2)) * sin(theta5)) * sin(theta6) - (-(-sin(theta1) * sin(theta3) + cos(theta1) * cos(theta2) * cos(theta3)) * sin(theta4) + sin(theta2) * cos(theta1) * cos(theta4)) * cos(theta6), 0, -sin(theta1) * sin(theta2), sin(theta1) * sin(theta2), -(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4), -(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4), -(((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * sin(theta6) + (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6), (((sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * cos(theta4) + sin(theta1) * sin(theta2) * sin(theta4)) * cos(theta5) - (sin(theta1) * sin(theta3) * cos(theta2) - cos(theta1) * cos(theta3)) * sin(theta5)) * sin(theta6) - (-(sin(theta1) * cos(theta2) * cos(theta3) + sin(theta3) * cos(theta1)) * sin(theta4) + sin(theta1) * sin(theta2) * cos(theta4)) * cos(theta6), 1, -cos(theta2), cos(theta2), sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4), sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4), -((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * sin(theta6) + (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6), ((-sin(theta2) * cos(theta3) * cos(theta4) + sin(theta4) * cos(theta2)) * cos(theta5) + sin(theta2) * sin(theta3) * sin(theta5)) * sin(theta6) - (sin(theta2) * sin(theta4) * cos(theta3) + cos(theta2) * cos(theta4)) * cos(theta6);

    return J;
}


// to ustvari kvaternion
Eigen::Quaterniond kvaternion(Eigen::MatrixXd M){

    Eigen::Matrix3d Nm;

    Nm << M(0,0),M(0,1),M(0,2),
          M(1,0),M(1,1),M(1,2),
          M(2,0),M(2,1),M(2,2);

    Eigen::Quaterniond kvaternion(Nm);

    return kvaternion;
}

Eigen::Vector3d RotMat2rpy(Eigen::MatrixXd M){

    Eigen::Matrix3d RotMat;
    RotMat << M(0, 0), M(0, 1), M(0, 2),
              M(1, 0), M(1, 1), M(1, 2),
              M(2, 0), M(2, 1), M(2, 2);

    Eigen::Quaterniond q(RotMat); 

    double roll = atan2(2*q.y()*q.w()-2*q.x()*q.z(), 1 - 2*q.y()*q.y()-2*q.z()*q.z());
    double pitch = atan2(2*q.x()*q.w() - 2*q.y()*q.z(), 1 - 2*q.x()*q.x() - 2*q.z()*q.z());
    double yaw = asin(2*q.x()*q.y() + 2*q.z()*q.w());

    Eigen::Vector3d rpy;
    rpy << pitch, roll, yaw;

    return rpy;
}
//koti ki so bližje joint limitom jih moras manj premikat 
//lahko nastavis sm kot konstante za en cikel?
double Gibljivost(Eigen::VectorXd Tzac, Eigen::VectorXd Tkon)
{
    //Eigen::VectorXd TrenutniKoti = *Koti;
    Eigen::VectorXd TrenutniKoti;

    for (int i = 0; i < 7; i++){
        TrenutniKoti(i)=seznam[i];
    }
    

    while(bool napaka = true){

        Eigen::VectorXd x = Tkon - Tzac;
        Eigen::VectorXd dx = x*0.1;

        Eigen::MatrixXd J = Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));
        Eigen::MatrixXd PsevdoJ = pseudoInverse(J);

        TrenutniKoti += PsevdoJ * dx;

        Eigen::MatrixXd Trenutno = DirektnaKinematika(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6));

        Eigen::VectorXd Orientacija = RotMat2rpy(Trenutno);

        Eigen::VectorXd Tt;

        Tt << Trenutno(0, 3), Trenutno(1, 3), Trenutno(2,3), Orientacija(0), Orientacija(1), Orientacija(2);

        Eigen::VectorXd razlika = Tkon - Tt;

        if (razlika.maxCoeff() <= 0.01){
            napaka = false;
            
            for(int i = 0; i<7; i++){
                seznam[i] = TrenutniKoti(i);//nastavi novo pozicijo kotov 
            }
        }
    }

    Eigen::MatrixXd VrednostPodKorenu = Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6)) * (Jakobi(TrenutniKoti(0), TrenutniKoti(1), TrenutniKoti(2), TrenutniKoti(3), TrenutniKoti(4), TrenutniKoti(5), TrenutniKoti(6)).transpose());

    double Yosihawa = sqrt(VrednostPodKorenu.determinant());
    return Yosihawa;

    //tu dodej da zapišeš novo vrednost kotov v seznam
}

int main()
{
    //nastavi pravilne zacetne vrednosti kotov 
    //Eigen::VectorXd zacetneVrednosti;
    //zacetneVrednosti << 0, 0, 0, 0, 0, 0, 0;
    //Eigen::VectorXd* Koti = &zacetneVrednosti;
    

    //std::cin.get();
    Eigen::VectorXd Tzac(6);
    Tzac << 0.088,0,0.926,0,0,0;
    Eigen::VectorXd Tkon(6);
    Tkon << 0.088, 0, 0.926, 0, 0, 0;
    std::cout<< DirektnaKinematika(0, 0, 0, 0, 0, 0, 0) << std::endl;
    std::cout << seznam[2] << std::endl;
    std::cout << Gibljivost(Tzac,Tkon);
    //std::cin.get();
}
