//
//  lab1_5.hpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//
/*
Реализовать алгоритм QR – разложения матриц в виде программы. На его основе разработать программу, реализующую QR – алгоритм решения полной проблемы собственных значений произвольных матриц, задавая в качестве входных данных матрицу и точность вычислений. С использованием разработанного программного обеспечения найти собственные значения матрицы.
*/

#ifndef lab1_5_hpp
#define lab1_5_hpp

#include <complex>
#include <utility>
#include "matrix.hpp"


double sign(double x) {
    if (x == 0) {
        return 0.0;
    }
    return x > 0 ? 1.0 : -1.0;
}

void QRDecompositionMethod(Matrix& A, Matrix& Q, Matrix& R) {
    Matrix H;
    Matrix E(A.getRow(), A.getColumn());
    E.unitary();
    // A = QR
    Q = E;
    R = A;
    for (int j = 0; j < A.getRow(); ++j) {
        std::vector<double> v(A.getRow(), 0.0);
        double norma = 0.0;
        
        v[j] = R[j][j];
        for (int i = 0; i < R.getRow(); ++i) {
            norma += R[i][j] * R[i][j];
        }
        norma = sqrt(norma);
        v[j] += sign(R[j][j]) * norma;
        for (int i = j + 1; i < R.getRow(); ++i) {
            v[i] = R[i][j];
        }

        double sum1 = 0.0;
        for (int i = 0; i < v.size(); ++i) {
            sum1 += v[i] * v[i];
        }
        Matrix sum2(v.size(), v.size());
        for (int i = 0; i < sum2.getRow(); ++i) {
            for (int j = 0; j < sum2.getColumn(); ++j) {
                sum2[i][j] = v[i] * v[j];
            }
        }
        H = E - ((2 / sum1) * sum2);
        Q = Q * H;
        R = H * R;
    }
}


std::pair<std::complex<double>, std::complex<double>> Roots(double a, double b, double c) {
    std::pair<std::complex<double>, std::complex<double>> res;
    std::complex<double> Discriminant = b * b - 4.0 * a * c;
    res.first = (-b + sqrt(Discriminant)) / (2 * a);
    res.second = (-b - sqrt(Discriminant)) / (2 * a);
    return res;
}


int32_t QRMethod(const Matrix& A, std::vector<std::complex<double>>& x, double accuracy) {
    Matrix Q, R, Ai = A;
    x.resize(Ai.getColumn());
    double criteria = 0.0;
    bool okay = true;
    int32_t iter = 0;

    while (okay) {
         QRDecompositionMethod(Ai, Q, R);
         Ai = R * Q;

         okay = false;
         for (int j = 0; j < Ai.getColumn(); ++j) {
             criteria = 0.0;
             for (int i = j + 1; i < Ai.getRow(); ++i) {
                 criteria += Ai[i][j] * Ai[i][j];
             }
             criteria = sqrt(criteria);

             if (criteria > accuracy) {
                 auto lambda = Roots(1.0, -Ai[j][j]-Ai[j+1][j+1], Ai[j][j]*Ai[j+1][j+1]-Ai[j][j+1]*Ai[j+1][j]);
                 auto x1 = abs(lambda.first - x[j]);
                 auto x2 = abs(lambda.second - x[j + 1]);
                 if (std::max(x1, x2) > accuracy) {
                     okay = true;
                 }
                 x[j] = lambda.first;
                 ++j;
                 x[j] = lambda.second;
             } else {
                 x[j] = Ai[j][j];
             }
         }
        ++iter;
     }
     return iter;
}

void InitMatrix5(Matrix& A) {
    A = Matrix(3, 3);
    A[0][0] = 9.0;
    A[0][1] = 0.0;
    A[0][2] = 2.0;

    A[1][0] = -6.0;
    A[1][1] = 4.0;
    A[1][2] = 4.0;

    A[2][0] = -2.0;
    A[2][1] = -7.0;
    A[2][2] = 5.0;
}

//void matrix_init(Matrix& A, int size){
//    A = Matrix(size, size);
//    for(int i = 0; i < size; ++i){
//        for(int j = 0; j < size; ++j){
//            std::cin >> A[i][j];
//        }
//    }
//}

bool lab1_5() {
    Matrix A;
    InitMatrix5(A);
    norma(A[0]);
    std::vector<std::complex<double>> x;
    double accuracy = 0.01;

    int iter = QRMethod(A, x, accuracy);
    
    std::cout << "S O L U T I O N:" << std::endl;
    for (int i = 0; i < x.size(); ++i) {
        std::cout << "l" << i + 1 << " = " << x[i].real() << " ";
        if (x[i].imag()) {
            if (x[i].imag() > 0) {
                std::cout << "+ ";
            }
            std::cout << x[i].imag() << "i";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Iterations: " << iter << std::endl;
    std::cout << iter << std::endl;
    return true;
}

#endif /* lab1_5_hpp */
