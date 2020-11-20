//
//  lab1_4.hpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#ifndef lab1_4_hpp
#define lab1_4_hpp

#include "matrix.hpp"

int32_t RotationMethod(Matrix& A, Matrix& Ui, std::vector<double>& x, double accuracy) {
    Matrix Ak = A;
    
    Ui = Matrix(A.getRow(), A.getColumn());
    Ui.unitary();
    
    x.resize(A.getColumn());
    
    Matrix U(A.getRow(), A.getColumn());
    
    int32_t iMax = 0, jMax = 0, iter = 0;
    double condition = 0.0, phi = 0.0, max;
    const double PI = 3.141592653589793;
    
    if (!A.isSimmetricalMatrix()){
        throw "Matrix not simmetric! Wrong!";
    }
    while(1) {
        max = 0.0;
        iMax = 0;
        jMax = 0;
        for (int i = 0; i < Ak.getRow(); ++i) {
            for (int j = i + 1; j < Ak.getColumn(); ++j) {
                if (max < abs(Ak[i][j])) {
                    max = abs(Ak[i][j]);
                    iMax = i;
                    jMax = j;
                }
            }
        }
        U.unitary();
        phi = Ak[iMax][iMax] != Ak[jMax][jMax] ? atan(2 * Ak[iMax][jMax] / (Ak[iMax][iMax] - Ak[jMax][jMax])) / 2 : PI/4;
        U[iMax][jMax] = -sin(phi);
        U[jMax][iMax] = sin(phi);
        U[iMax][iMax] = cos(phi);
        U[jMax][jMax] = cos(phi);

        Ui = Ui * U;

        U.transpose();
        Ak = U * Ak;
        U.transpose();

        Ak = Ak * U;
        
        condition = 0;
        for (int i = 0; i < Ak.getRow(); ++i) {
            for (int j = i + 1; j < Ak.getColumn(); ++j) {
                condition += Ak[i][j] * Ak[i][j];
            }
        }
        condition = sqrt(condition);
        if (condition < accuracy) {
            break;
        }
        ++iter;
    }
    for (int i = 0; i < Ak.getRow(); ++i){
        x[i] = Ak[i][i];
    }
    return iter;
}

void InitMatrix4(Matrix& A) {
    A = Matrix(3, 3);
    A[0][0] = -6.0;
    A[0][1] = 6.0;
    A[0][2] = -8.0;

    A[1][0] = 6.0;
    A[1][1] = -4.0;
    A[1][2] = 9.0;

    A[2][0] = -8.0;
    A[2][1] = 9.0;
    A[2][2] = -2.0;
}

bool lab1_4() {
    Matrix A, U;
    InitMatrix4(A);
    std::cout << A << std::endl;
    std::vector<double> x;
    double accuracy = 0.001;
    int32_t iter = RotationMethod(A, U, x, accuracy);
    
    std::cout << "Sobstv values:" << std::endl;
    for (int i = 0; i < x.size(); ++i){
        std::cout << "x" << i + 1 << " = " << x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Iterations = " << iter << std::endl;
    std::cout << "Sobstv vectors:" << std::endl;
    std::vector<double> vec_x(U.getColumn());
    std::vector<double> solu2(U.getColumn());
    for (int j = 0; j < U.getColumn(); ++j){
              std::cout << "x" << j + 1 << " = " << std::endl;
              for (int i = 0; i < U.getRow(); ++i){
                  std::cout << U[i][j] << std::endl;
              }
           std::cout << std::endl;

    }
    
    return true;
}

#endif /* lab1_4_hpp */
