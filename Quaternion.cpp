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
Eigen::Quaterniond kvaternion;
Eigen::Vector4d vektor;

int main(){

    kvaternion = Eigen::Quaterniond(3, 5, 6, 7);
    vektor << kvaternion.w(),kvaternion.x(),kvaternion.y(),kvaternion.z();
    cout << vektor;
}