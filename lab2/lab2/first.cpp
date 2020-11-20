//
//  main.cpp
//  lab2
//
//  Created by Лина Вельтман on 06/03/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//
/*
 2.1. Реализовать методы простой итерации и Ньютона решения нелинейных уравнений в виде программ, задавая в качестве входных данных точность вычислений. С использованием разработанного программного обеспечения найти положительный корень нелинейного уравнения (начальное приближение определить графически). Проанализировать зависимость погрешности вычислений от количества итераций.

 */

#include <iostream>
#include <cmath>
#include <vector>


const double epsilon = 0.001;

double foo(double x) {
    return (pow(2, x) + x * x - 2);
//    return pow(e, 2*x) + 3*x-4;
}


double firstDiff(double x) {
//    return 2*pow(e,2*x)+3;
    return (pow(2, x) * log(2) + 2 * x);
}


double secondDiff(double x) {
//    return 4 * pow(e, 2*x);
    return (pow(2, x) * log(2) * log(2) + 2);
}


int NewtonsMethod(double a, double b, double eps, double& ans) {
    int k = 0;
    double xk_1 = /*0.6*/a + eps;
    double xk = xk_1;
    if (foo(a) * foo(b) < 0) {
        double criteria = eps;
        while (criteria >= eps) {
            if (foo(xk_1) * secondDiff(xk_1) > 0) {
                xk = xk_1;
                xk_1 = xk - (foo(xk)) / firstDiff(xk);
                ++k;
                criteria = abs(xk_1 - xk);
            } else {
                if (xk_1 < b) {
                    xk_1 += eps;
                } else {
                    throw "Try to find another length [a;b]";
                }
            }
        }
    } else {
        throw "Try to find another length [a;b]";
    }
    ans = xk_1;
    return k + 1;
}

double phi(double x) {
    // + -
    return sqrt(2 - pow(2, x));
//    return log(4-3*x)/2;
}


double phiDiff(double x) {
    return (-pow(2, x)* log(2)) / (2 * sqrt(2 - pow(2, x)));
//    return -3/(2*(4-3*x));
}

// dihotomy methood for search supr:
double get_sup(double a, double b) {
    double F1, F2;
    while (abs(b - a) >= epsilon) {
        double x = (a + b) / 2.0;
        F2 = phiDiff(x + epsilon);
        F1 = phiDiff(x - epsilon);
        if (F1 < F2) {
            a = x;
        } else {
            b = x;
        }
    }
    return phiDiff((b + a) / 2.0);
}

double newton_method(double x0, double alfa, int& itter){
    double x_j, x_k = x0;
    itter = 0;
    do{
        x_j = x_k;
        x_k -= foo(x_j) / firstDiff(x_j);
        ++itter;
    } while(abs(x_k - x_j) >= alfa);
    return x_k;
}

double itteration_method(double x0, double alfa, int& itter, double a, double b) {
    double x_j, x_k = x0;
    // search q here:
    double q = get_sup(a, b);
    q /= (1 - q);

    itter = 0;
    do {
        x_j = x_k;
        x_k = phi(x_j);
        ++itter;
    } while(q * abs(x_k - x_j) >= alfa);
    return x_k;
}

int main(int argc, const char * argv[]) {
    double accuracy, x0;
    double a = 0.0, b = 1.0, ans1 = 0.0, ans2 = 0.0;
    //double x0 = (a + b) / 2;
    std::cin >> accuracy;
    std::cin >> x0;
    int iterNewton, iterSimple;
    
    ans1 = newton_method(x0, accuracy, iterNewton);
    ans2 = itteration_method(x0, accuracy, iterSimple, a, b);
    std::cout << "ANSWER:" << std::endl;
    std::cout << "Newton's method" << std::endl;
    std::cout << "x = " << ans1 << std::endl;
    std::cout << "Iterations: " << iterNewton << std::endl;
    std::cout << "Simple iterations method" << std::endl;
    std::cout << "x = " << ans2 << std::endl;
    std::cout << "Iterations: " << iterSimple << std::endl;
    return 0;
}
