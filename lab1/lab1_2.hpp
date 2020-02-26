//
//  main.cpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. mll rights reserved.
//
#ifndef lab1_2_hpp
#define lab1_2_hpp

#include "matrix.hpp"

/*
     Метод прогонки (англ. tridiagonal matrix algorithm)
     или алгоритм Томаса (англ. Thomas algorithm)
     используется для решения систем линейных уравнений вида
     Ax = F, где A — трёхдиагональная матрица.
 */

void ThomasAlgorithm(std::vector<std::vector<double>>& koef, std::vector<double>& d, std::vector<double>& x) {
    std::vector<double> p(d.size()), q(d.size());
    p[0] = -koef[0][2] / koef[0][1];
    q[0] = d[0] / koef[0][1];
    int n = (int) x.size();
    for (int i = 1; i < n; ++i) {
        p[i] = -koef[i][2] / (koef[i][1] + koef[i][0] * p[i - 1]);
        q[i] = (d[i] - koef[i][0] * q[i - 1]) / (koef[i][1] + koef[i][0] * p[i - 1]);
    }
    x[n - 1] = q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = p[i] * x[i + 1] + q[i];
    }
}

void InitMatrix(Matrix& m) {
    m = Matrix(5, 5);
    m[0][0] = 15.0;
    m[0][1] = 8.0;
    m[0][2] = 0.0;
    m[0][3] = 0.0;
    m[0][4] = 0.0;

    m[1][0] = 2.0;
    m[1][1] = -15.0;
    m[1][2] = 4.0;
    m[1][3] = 0.0;
    m[1][4] = 0.0;

    m[2][0] = 0.0;
    m[2][1] = 4.0;
    m[2][2] = 11.0;
    m[2][3] = 5.0;
    m[2][4] = 0.0;

    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = -3.0;
    m[3][3] = 16.0;
    m[3][4] = -7.0;

    m[4][0] = 0.0;
    m[4][1] = 0.0;
    m[4][2] = 0.0;
    m[4][3] = 3.0;
    m[4][4] = 8.0;
}

void InitVectorB(std::vector<double>& b) {
    b[0] = 92.0;
    b[1] = -84.0;
    b[2] = -77.0;
    b[3] = 15.0;
    b[4] = -11.0;
}

bool lab1_2() {
    Matrix m;
    InitMatrix(m);
    std::vector<double> b(m.getRow()), x(m.getRow(), 0.0);
    InitVectorB(b);
    
    try {
        if (!m.isTridiagonalMatrix()) {
            throw 3;
        }
        std::vector<std::vector<double>> tmp;
        m.getKoefficients(tmp);
        ThomasAlgorithm(tmp, b, x);
        for (int i = 0; i < x.size(); ++i) {
            std::cout << "X" << i + 1 << " = " << x[i] << " ";
        }
        std::cout << "\n";
        
    } catch (int32_t i) {
        std::cout << "Error №" << i << ": matrix must be tridiagonal!" << std::endl;
        return false;
    }
    
    return true;
}

#endif /* lab1_2_hpp */
