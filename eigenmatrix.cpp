#include <iostream>
#include <Eigen/Dense>
#include <complex>

int main(){
    int n = 2; // number of nodes
    Eigen::MatrixXcf matrixC(n,n);
    std::cout << "creation of a matrix with " << n << " nodes" << std::endl;

    std::cout << std::endl << "set all values to 0: " << std::endl;
    matrixC.setZero();
    std::cout << matrixC << std::endl;


    std::cout << std::endl << "change required values: " << std::endl;

    const std::complex<float> Re(1.0f, 0.0f); // set real part
    const std::complex<float> Im(0.0f, 1.0f); // set imaginary part

    std::complex<float> a, b;
    a = 3;
    b = 1;

    matrixC(0,0) = a * Re + b * Im;
    matrixC(0,1) = 3.0f + 1.0f * Im;
    matrixC(1,0) = 4.0f + 0.0f * Im;
    matrixC(1,1) = 1.0f + 2.0f * Im;

    std::cout << matrixC << std::endl;


    Eigen::MatrixXcf Cinv;
    Cinv = matrixC.inverse();
    std::cout << std::endl << "the inverse of the matrix is: " << std::endl;
    std::cout << Cinv << std::endl;

    std::cout << std::endl << "check the accuracy by multiplying matrix by its inverse: " << std::endl;
    std::cout << matrixC * Cinv << std::endl;

    std::cout << std::endl << "column vector: " << std::endl;
    Eigen::MatrixXcf columnVector(n,1);
    columnVector << 1, 1;
    std::cout << columnVector << std::endl;

    std::cout << std::endl << "multiplication of the inverse matrix with a column vector: " << std::endl;
    Eigen::MatrixXcf solution(n,1);
    solution = Cinv * columnVector;
    std::cout << solution << std::endl << std::endl;

    for(int i = 0; i < n; i++){
        std::cout << "solution at (" << i << ",0): " << solution(i,0) << std::endl;
    }
}
