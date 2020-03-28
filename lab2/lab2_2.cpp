//
//  lab2_2.cpp
//  lab2
//
//  Created by Лина Вельтман on 10/03/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <limits>
//#include "matrix.hpp"



std::pair<double, double> FuncSystem(double x1, double x2) {
    return std::make_pair(
        //x1 * x1 + x2 * x2 - 4,
        0.1 * x1 * x1 + x1 + 0.2 * x2 * x2 - 0.3,
        0.2 * x2 * x2 + x2 - 0.1 * x1 * x2 - 0.7
        //x1 - exp(x2) + 2
    );
}

double FuncSystemFi1(double x1) {
    //return sqrt(- x1 * x1 + 4);
    return sqrt((0.3-0.1*x1*x1-x1)/0.2);
}

double FuncSystemFi2(double x2) {
//    return exp(x2) - 2;
    return ((0.2*x2*x2+x2-0.7) / (0.1*x2));
}

double FuncSystemFi1dx(double x1) {
//    return - x1 / sqrt(- x1 * x1 + 4);
    return (-0.5*x1 - 2.5)/(sqrt(1.5-0.5*x1*x1-5*x1));
}

double FuncSystemFi2dx(double x2) {
    //return exp(x2);
    return ((10*(0.4*x2+1))/x2) - ((10*(0.2*x2*x2+x2-0.7))/(x2*x2));
}

std::pair<double, double> FiSystem(double x1, double x2) {
    return std::make_pair(
        FuncSystemFi1(x1),
        FuncSystemFi2(x2)
    );
}

std::pair<double, double> FiSystemDx(double x1, double x2) {
    return std::make_pair(
        FuncSystemFi1dx(x1),
        FuncSystemFi2dx(x2)
    );
}

double sup(double a, double b, double eps) {
    std::pair<double, double> F1, F2;
    while (abs(b - a) >= eps){
        double x = (a + b) / 2.0;
        F1 = FiSystemDx(x - eps, x - eps);
        F2 = FiSystemDx(x + eps, x + eps);
        if (F1.first < F2.second && F1.second < F2.second) {
            a = x;
        } else {
            b = x;
        }
    }
    F1 = FiSystemDx((b + a) / 2.0, (b + a) / 2.0);
    return std::max(F1.first, F1.second);
}

std::pair<double, double> NewtonMethod(double x1, double x2, double eps) {
    Matrix j(2, 2);
    Matrix jr(2, 2);
    double jdet;
    std::pair<double, double> func_curr = FuncSystem(x1, x2);

    j[0][0] = (FuncSystem(x1 + eps, x2).first - func_curr.first) / eps;
    j[0][1] = (FuncSystem(x1, x2 + eps).first - func_curr.first) / eps;
    j[1][0] = (FuncSystem(x1 + eps, x2).second - func_curr.second) / eps;
    j[1][1] = (FuncSystem(x1, x2 + eps).second - func_curr.second) / eps;
    jdet = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    jr[0][0] = j[1][1] / jdet;
    jr[0][1] = -j[0][1] / jdet;
    jr[1][0] = -j[1][0] / jdet;
    jr[1][1] = j[0][0] / jdet;

    return std::make_pair(
        x1 - jr[0][0] * func_curr.first - jr[0][1] * func_curr.second,
        x2 - jr[1][0] * func_curr.first - jr[1][1] * func_curr.second
    );
}

int doAlgo(std::pair<double, double>& xk, std::pair<double, double>& xk1, int method) {

    const double eps = 0.0001;/*(double) std::numeric_limits<double>::epsilon() * 10000000.;*/
    xk = std::make_pair(0.25, 0.75);
    xk1 = std::make_pair(0.25, 0.75);
    int iter = 0;
    double q = sup(xk.first, xk.second, eps);
    if (method) {
        do {
            ++iter;
            xk = xk1;
            xk1 = NewtonMethod(xk.first, xk.second, eps);
        } while (std::max(abs(xk1.first - xk.first), abs(xk1.second - xk.second)) >= eps);
    } else {
        do {
            ++iter;
            xk = xk1;
            xk1 = FiSystem(xk.first, xk.second);
        } while (q/(1 - q) * std::max(abs(xk1.first - xk.first), abs(xk1.second - xk.second)) >= eps);
    }
    return iter;
}

int main() {
    std::pair<double, double> xk, xk1;
    std::cout << "ANSWER:" << std::endl;
    
    std::cout << "Newton's method" << std::endl;
    int iterNewton = doAlgo(xk, xk1, 1);
    std::cout << "X = (" << xk1.first << ", " << xk1.second << ")" << std::endl;

//    std::cout << "x = (" << x_n[0];
//    for(int i = 1; i < n; ++i){
//        std::cout << ", " << x_n[i];
//    }
//    std::cout << ")" << std::endl;
    std::cout << "Iterations: " << iterNewton << std::endl;
    
    std::cout << "Simple iterations method" << std::endl;
    int iterSimple = doAlgo(xk, xk1, 0);
    std::cout << "X = (" << xk1.first << ", " << xk1.second << ")" << std::endl;

//    std::cout << "x = (" << x_n[0];
//    for(int i = 1; i < n; ++i){
//        std::cout << ", " << x_i[i];
//    }
//    std::cout << ")" << std::endl;
    std::cout << "Iterations: " << iterSimple << std::endl;
    return 0;
}
