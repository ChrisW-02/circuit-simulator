#include <iostream>
#include <complex>
#include <string>
#include <math.h>

#include "Matrix.h"

// a simple function to print matrix
template <class T>
void PrintMatrix(Matrix<T> matrix){
    int nRows = matrix.GetNumRows();
    int nCols = matrix.GetNumCols();
    for(int row = 0; row < nRows; ++row){
        for(int col = 0; col < nCols; ++col){
            std::cout << std::setprecision(3) << matrix.GetElement(row, col) << "\t";
        }
        std::cout << std::endl;
    }
}

int main(){
    float matData[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix<float> testingMatrix(3, 3, matData);

    std::cout << std::endl << "----------------------------------------" << std::endl;
    std::cout << "3x3 matrix test (testingMatrix)" << std::endl;
    PrintMatrix(testingMatrix);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test matrix inversion of a square" << std::endl;
    float invertTestData[9] = {2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};
    Matrix<float> invertTest(3, 3, invertTestData);
    Matrix<float> invertResult = invertTest;
    invertResult.Inverse();
    std::cout << "From:" << std::endl;
    PrintMatrix(invertTest);
    std::cout << "To:" << std::endl;
    PrintMatrix(invertResult);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test multiplication of a matrix by its inverse" << std::endl;
    std::cout << "using invertTest * invertResult:" << std::endl;
    Matrix<float> invertAccuracy = invertTest * invertResult;
    PrintMatrix(invertAccuracy);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test multiplication of a column vector and an inverse matrix" << std::endl;
    float columnData2[3] = {1, 2, 3};
    Matrix<float> testColumn2(3, 1, columnData2);
    std::cout << "Column vector: " << std::endl;
    PrintMatrix(testColumn2);
    std::cout << "Square matrix: " << std::endl;
    PrintMatrix(invertResult);
    std::cout << "Square matrix * Column vector = " << std::endl;
    PrintMatrix(invertResult * testColumn2);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Testing complete" << std::endl;
    std::cout << "------------------------------" << std::endl;

    // test case examples below

    /*
    // create an instance of the Matrix class
    // this will contain a simple 2D 3x3 matrix
    float simpleData[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix<float> testMatrix(3, 3, simpleData);

    // extract and print the elements of testMatrix
    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "3x3 matrix test (testMatrix)" << std::endl;
    PrintMatrix(testMatrix);

    // test element retrieval
    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test Element Retrieval" << std::endl;
    std::cout << "Element (0,0) = " << testMatrix.GetElement(0,0) << std::endl;
    std::cout << "Element (1,0) = " << testMatrix.GetElement(1,0) << std::endl;
    std::cout << "Element (2,0) = " << testMatrix.GetElement(2,0) << std::endl;
    std::cout << "Element (1,1) = " << testMatrix.GetElement(1,1) << std::endl;
    std::cout << "Element (2,2) = " << testMatrix.GetElement(2,2) << std::endl;
    std::cout << "Element (3,3) = " << testMatrix.GetElement(3,3) << std::endl;
    std::cout << "Element (5,5) = " << testMatrix.GetElement(5,5) << std::endl;

    // test matrix multiplication
    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test Matrix Multiplication of column vector and matrix" << std::endl;
    double columnData[3] = {1, 1, 1};
    double squareData[9] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    Matrix<double> testColumn(3, 1, columnData);
    Matrix<double> squareMatrix(3, 3, squareData);
    std::cout << "Column vector: " << std::endl;
    PrintMatrix(testColumn);
    std::cout << "Square matrix: " << std::endl;
    PrintMatrix(squareMatrix);
    std::cout << "Column vector * Square matrix = " << std::endl;
    PrintMatrix(testColumn * squareMatrix);
    std::cout << "Square matrix * Column vector = " << std::endl;
    PrintMatrix(squareMatrix * testColumn);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test formation of identity matrix" << std::endl;
    Matrix<double> identityTest(5,5);
    identityTest.SetToIdentity();
    PrintMatrix(identityTest);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test joining of two matrices" << std::endl;
    Matrix<double> bigSquare(5,5);
    bigSquare.Join(identityTest);
    PrintMatrix(bigSquare);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test matrix inversion of a square" << std::endl;
    double invertTestData[9] = {2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};
    Matrix<double> invertTest(3, 3, invertTestData);
    Matrix<double> invertResult = invertTest;
    invertResult.Inverse();
    std::cout << "From:" << std::endl;
    PrintMatrix(invertTest);
    std::cout << "To:" << std::endl;
    PrintMatrix(invertResult);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test multiplication of a matrix by its inverse" << std::endl;
    std::cout << "using invertTest * invertResult:" << std::endl;
    Matrix<double> invertAccuracy = invertTest * invertResult;
    PrintMatrix(invertAccuracy);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test multiplication of a column vector and an inverse matrix" << std::endl;
    double columnData2[3] = {1, 2, 3};
    Matrix<double> testColumn2(3, 1, columnData2);
    std::cout << "Column vector: " << std::endl;
    PrintMatrix(testColumn2);
    std::cout << "Square matrix: " << std::endl;
    PrintMatrix(invertResult);
    std::cout << "Square matrix * Column vector = " << std::endl;
    PrintMatrix(invertResult * testColumn2);

    std::cout << std::endl << "------------------------------" << std::endl;
    std::cout << "Test bigger matrix inversion" << std::endl;
    double invertTestData2[25] =
			{2.0, 3.0, 4.0, 5.0, 6.0,
			 1.0, 2.0, 3.0, 4.0, 5.0,
			 9.0, 5.0, 3.0, 2.0, 6.0,
			 2.0, 4.0, 6.0, 5.0, 1.0,
			 1.0, 7.0, 5.0, 2.0, 3.0};
    Matrix<double> invertTest2(5, 5, invertTestData2);
    Matrix<double> invertResult2 = invertTest2;
    invertResult2.Inverse();
    std::cout << "From:" << std::endl;
    PrintMatrix(invertTest2);
    std::cout << "To:" << std::endl;
    PrintMatrix(invertResult2); */
}
