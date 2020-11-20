//
//  lab1_3.hpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#ifndef lab1_3_hpp
#define lab1_3_hpp

#include "matrix.hpp"

double norma(const std::vector<double>& m) {
    double res = 0.0;
    for (int i = 0; i < m.size(); ++i) {
        res += m[i] * m[i];
    }
    return sqrt(res);
}

std::vector<double> Substraction(std::vector<double>& x, std::vector<double>& new_x) {
    std::vector<double> vec = x;
    for (int i = 0; i < x.size(); ++i) {
        vec[i] -= new_x[i];
    }
    return vec;
}

std::vector<double> Addition(const std::vector<double>& first, const std::vector<double>& second) {
    int32_t size = (int32_t) std::max(first.size(), second.size());
    std::vector<double> res(size, 0.0);
    for (int i = 0; i < size; ++i) {
        res[i] = first[i] + second[i];
    }
    return res;
}

double Checker(std::vector<double>& x, std::vector<double>& old_x) {
    std::vector<double> vec = Substraction(x, old_x);
    return norma(vec);
}

void ToEquivalentForm(const Matrix& A, Matrix& alpha, std::vector<double>& beta) {
    if (!A.isSquareMatrix()) {
       throw "It's not a square matrix!\n";
    }
    for (int i = 0; i < alpha.getRow(); ++i) {
        if (!A[i][i]) {
            throw "Main diagonal consists 0! Swap the appropriate equation with any other equation!\n";
        }
        for (int j = 0; j < alpha.getRow(); ++j) {
            alpha[i][j] = i != j ? -A[i][j] / A[i][i] : 0.0;
        }
        beta[i] /= A[i][i];
    }
}

int32_t SimpleIterationsMethod(const Matrix A, const std::vector<double> b, std::vector<double>& x, double accuracy) {
    Matrix alpha = A;
    x.resize(b.size());
    std::vector<double> prev_x(b.size(), 0.0),  beta = b;
    
    double coeff = 0.0;
    
    ToEquivalentForm(A, alpha, beta);
    
    x = beta;
    
    coeff = alpha.get_norm();
    if (coeff < 1.) {
        coeff /= (1 - coeff);
    } else {
        coeff = 1.0;
    }
    
    int32_t iter = 0;
    std::cout << coeff * Checker(x, prev_x) << std::endl;
    for (iter = 0; coeff * Checker(x, prev_x) > accuracy; ++iter) {
        prev_x = x;
        x = Addition(beta, alpha * prev_x);
    }
    return iter;
}

int32_t ZeidelMethod(const Matrix A, const std::vector<double> b, std::vector<double>& x, double accuracy) {
    Matrix alpha = A;
    x.resize(b.size());
    std::vector<double> prev_x(b.size(), 0.0), beta = b;
    
    double coeff = 0.0;
    
    ToEquivalentForm(A, alpha, beta);
    
    x = beta;
    
    coeff = alpha.get_norm();
    if (coeff < 1.0) {
        coeff = alpha.get_upper_norma() / (1 - coeff);
    } else {
        coeff = 1.0;
    }
    
    int32_t iter = 0;
    for (iter = 0; coeff * Checker(x, prev_x) > accuracy; ++iter) {
         x.swap(prev_x);
         x = beta;
         for (int i = 0; i < alpha.getRow(); ++i) {
             for (int j = 0; j < i; ++j){
                 x[i] += x[j] * alpha[i][j];
             }
             for (int j = i; j < alpha.getColumn(); ++j) {
                 x[i] += prev_x[j] * alpha[i][j];
             }
         }
     }
    return iter;
}


void InitMatrix3(Matrix& m) {
    m = Matrix(4, 4);
    m[0][0] = 29.0;
    m[0][1] = 8.0;
    m[0][2] = 9.0;
    m[0][3] = -9.0;

    m[1][0] = -7.0;
    m[1][1] = -25.0;
    m[1][2] = 0.0;
    m[1][3] = 9.0;

    m[2][0] = 1.0;
    m[2][1] = 6.0;
    m[2][2] = 16.0;
    m[2][3] = -2.0;

    m[3][0] = -7.0;
    m[3][1] = 4.0;
    m[3][2] = -2.0;
    m[3][3] = 17.0;
}

void InitVectorB3(std::vector<double>& b) {
    b.resize(4);
    b[0] = 197.0;
    b[1] = -226.0;
    b[2] = -95.0;
    b[3] = -58.0;
}

void Answer(std::vector<double>& x, int32_t iter) {
    std::cout << "S O L U T I O N:" << std::endl;
    for (int i = 0; i < x.size(); ++i) {
        std::cout << "x" << i + 1 << " = " << x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Iterations: " << iter << std::endl;
}

bool lab1_3() {
    double accuracy = 0.01;
    
    Matrix m;
    InitMatrix3(m);
    Matrix alpha = m;
    std::cout << m << std::endl;
    std::vector<double> b(m.getRow()), x;
    InitVectorB3(b);
    
    int iter = SimpleIterationsMethod(alpha, b, x, accuracy);
    std::cout << "SimpleIterationsMethod" << std::endl;
    Answer(x, iter);

    std::cout << "ZeidelMethod" << std::endl;
    iter = ZeidelMethod(alpha, b, x, accuracy);
    Answer(x, iter);

    return true;
}

#endif /* lab1_3_hpp */
