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

using namespace std;

int main(){
    Eigen::Vector3d e1;
    Eigen::Vector3d e2;

    e1 << 1,2,3;
    e2 << 4,5,6;

    cout << e1 << endl;
    cout << "----" << endl;
    cout << e2 - e1 <<endl;
    cout << "----" << endl;
    cout << e1.cross(e2) << endl;
}