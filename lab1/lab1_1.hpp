
//
//  lab1_1.h
//  lab1
//
//  Created by Лина Вельтман on 04/03/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#ifndef lab1_1_h
#define lab1_1_h


#include "matrix.hpp"

void NewSolutionUsingLUDecomposition(Matrix& L, Matrix& U, std::vector<double>& b, std::vector<double>& x) {
    x.resize(U.getColumn());
    std::vector<double> z(L.getColumn());
    z[0] = b[0];
    for(int i = 1; i < L.getColumn(); ++i) {
        z[i] = b[i];
        for(int j = 0; j < i; ++j){
            z[i] -= (L[i][j] * z[j]);
        }
    }
    
    for (int i = U.getColumn() - 1; i >= 0; --i) {
        x[i] = z[i];
        for(int j = i + 1; j < U.getColumn(); ++j) {
            x[i] -= (x[j] * U[i][j]);
        }
        x[i] /= U[i][i];
    }
}

void Gauss(Matrix& L, Matrix& U, Matrix& matrix, std::vector<double>& b, std::vector<double>& x) {
   
    U = matrix;
    L = Matrix(matrix.getRow(), matrix.getColumn());
    for (int k = 0; k < U.getRow(); ++k) {
        if (U[k][k] == 0.0){
            throw "LU not exist! Try simple gauss method or enter lines in another order!";
        }
        L[k][k] = 1.0;
        for (int i = k + 1; i < U.getRow(); ++i) {
            double m_ik = U[i][k] / U[k][k];
            L[i][k] = m_ik;

            for (int j = k; j < U.getRow(); ++j) {
                U[i][j] = U[i][j] - m_ik * U[k][j];
            }
        }
    }
    NewSolutionUsingLUDecomposition(L, U, b, x);
}
    

void InvertMartix(Matrix& L, Matrix& U, Matrix& invA, std::vector<double>& b) {
    std::vector<double> tmp(b.size());
    for (int i = 0; i < invA.getRow(); ++i) {
        tmp[i] = 1.0;
        NewSolutionUsingLUDecomposition(L, U, tmp, invA[i]);
        tmp[i] = 0.0;
    }
    invA.transpose();
}

void InitMatrix1(Matrix& A) {
    A = Matrix(4, 4);
    A[0][0] = 1.0;
    A[0][1] = -5.0;
    A[0][2] = -7.0;
    A[0][3] = 1.0;

    A[1][0] = 1.0;
    A[1][1] = -3.0;
    A[1][2] = -9.0;
    A[1][3] = -4.0;

    A[2][0] = -2.0;
    A[2][1] = 4.0;
    A[2][2] = 2.0;
    A[2][3] = 1.0;

    A[3][0] = -9.0;
    A[3][1] = 9.0;
    A[3][2] = 5.0;
    A[3][3] = 3.0;
}

void InitVectorB1(std::vector<double>& b) {
    b.resize(4);
    b[0] = -75.0;
    b[1] = -41.0;
    b[2] = 18.0;
    b[3] = 29.0;
}

void Solution(Matrix& L, Matrix& U, Matrix& invA, std::vector<double>& x, double determinant) {
    std::cout << "S O L U T I O N:" << std::endl;
    std::cout << L << std::endl;
    std::cout << U << std::endl;
    std::cout << "Inverse matrix:" << std::endl;
    std::cout << invA << std::endl;
    std::cout << "x solution:" << std::endl;
    for (int i = 0; i < x.size(); ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Determinant: " << determinant << std::endl;
}

bool lab1_1() {
    Matrix L, U, A;
    InitMatrix1(A);
    Matrix invA(A.getColumn(), A.getRow());

    std::vector<double> b, x;
    InitVectorB1(b);
    
    Gauss(L, U, A, b, x);
    InvertMartix(L, U, invA, b);
    double determinant = 1.0;

    for (int i = 0; i < U.getRow(); ++i) {
         determinant *= U[i][i];
     }

    Solution(L, U, invA, x, determinant);

    return true;
}


#endif /* lab1_1_h */
