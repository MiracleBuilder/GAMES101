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
    /*
        给定一个点 P=(2,1), 将该点绕原点先逆时针旋转 45度，再平移 (1,2), 
        计算出变换后点的坐标（要求用齐次坐标进行计算）。
    */
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