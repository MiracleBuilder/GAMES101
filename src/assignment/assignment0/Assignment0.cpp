#define _USE_MATH_DEFINES
#include <cmath>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>

#define EHD 180.0;
#define PI 3.1415926f;

using namespace std;

int main()
{
    // counterclockwise 45 degree
    float alpha = M_PI / 4;

    float cosAlpha = std::cos(alpha);
    float sinAlpha = std::sin(alpha);
    int tx = 1, ty = 2;
    Eigen::Matrix3f rotateMatrix;
    rotateMatrix << cosAlpha, -sinAlpha, (float)tx,
        sinAlpha, cosAlpha, (float)ty,
        0, 0, 1;

    std::cout << rotateMatrix << std::endl;
    Eigen::Vector3f position(2.0f, 1.0f, 1.0f);
    position = rotateMatrix * position;
    std::cout << "after rotation: " << std::endl
              << position << std::endl;

    Eigen::Vector2f v;
    v << 2.0f, 1.0f;

    Eigen::Rotation2Df t(alpha);
    t.toRotationMatrix();
    printf("\nUsing an Affine2f\n");
    std::cout << "\nthe rotation Matrix is \n"
              << t.matrix() << std::endl;
              
    Eigen::Vector2f rotatedVect = t * v;
    std::cout << "the rotated vector is \n"
              << rotatedVect << std::endl;
}