////
////  main.cpp
////  lab2
////
////  Created by Лина Вельтман on 06/03/2020.
////  Copyright © 2020 Lina Veltman. All rights reserved.
////
///*
// 2.1. Реализовать методы простой итерации и Ньютона решения нелинейных уравнений в виде программ, задавая в качестве входных данных точность вычислений. С использованием разработанного программного обеспечения найти положительный корень нелинейного уравнения (начальное приближение определить графически). Проанализировать зависимость погрешности вычислений от количества итераций.
//
// */
//
//#include <iostream>
//#include <cmath>
//#include <vector>
//
////const double e = 2.71828182846;
//double foo(double x) {
//    return (pow(2, x) + x * x - 2);
////    return pow(e, 2*x) + 3*x-4;
//}
//
//
//double firstDiff(double x) {
////    return 2*pow(e,2*x)+3;
//    return (pow(2, x) * log(2) + 2 * x);
//}
//
//
//double secondDiff(double x) {
////    return 4 * pow(e, 2*x);
//    return (pow(2, x) * log(2) * log(2) + 2);
//}
//
//
//int NewtonsMethod(double a, double b, double eps, double& ans) {
//    int k = 0;
//    double xk_1 = /*0.6*/a + eps;
//    double xk = xk_1;
//    if (foo(a) * foo(b) < 0) {
//        double criteria = eps;
//        while (criteria >= eps) {
//            if (foo(xk_1) * secondDiff(xk_1) > 0) {
//                xk = xk_1;
//                xk_1 = xk - (foo(xk)) / firstDiff(xk);
//                ++k;
//                criteria = abs(xk_1 - xk);
//            } else {
//                if (xk_1 < b) {
//                    xk_1 += eps;
//                } else {
//                    throw "Try to find another length [a;b]";
//                }
//            }
//        }
//    } else {
//        throw "Try to find another length [a;b]";
//    }
//    ans = xk_1;
//    return k + 1;
//}
//
//double phi(double x) {
//    // + -
//    return sqrt(2 - pow(2, x));
////    return log(4-3*x)/2;
//}
//
//
//double phiDiff(double x) {
//    return (-pow(2, x)* log(2)) / (2 * sqrt(2 - pow(2, x)));
////    return -3/(2*(4-3*x));
//}
//
//
//int LyambdaSolution(double xk, double& xk_1, double eps) {
//    int k = 0;
//    double lyambda = eps;
//    while (true)
//        if (abs(1 - lyambda * firstDiff(xk_1)) < 1) {
//            double criteria = eps;
//            while (criteria >= eps) {
//                xk = xk_1;
//                xk_1 = phi(xk);
//                ++k;
//                criteria = abs(xk_1 - xk);
//            }
//            break;
//        } else {
//            lyambda += eps;
//        }
//    return k + 1;
//}
//
//int SimpleIterationsMethod(double a, double b, double eps, double& ans) {
//    double x0 = (a + b) / 2, xk;
//    double xk_1 = xk = x0;
//    int k = 0;
//    if (phi(x0) >= a && phi(x0) <= b) {
//        if (abs(phiDiff(x0)) < 1) {
//            double criteria = eps;
//            while (criteria >= eps) {
//                xk = xk_1;
//                xk_1 = phi(xk);
//                ++k;
//                criteria = abs(xk_1 - xk);
//            }
//        } else {
//            k = LyambdaSolution(xk, xk_1, eps);
//        }
//    } else {
//        k = LyambdaSolution(xk, xk_1, eps);
//    }
//    ans = xk_1;
//    return k + 1;
//}
//int main(int argc, const char * argv[]) {
//    double accuracy = 0.001;
//    double a = 0.0, b = 1.0, ans1 = 0.0, ans2 = 0.0;
//    int iterNewton = NewtonsMethod(a, b, accuracy, ans1);
//    int iterSimple = SimpleIterationsMethod(a, b, accuracy, ans2);
//
//    std::cout << "ANSWER:" << std::endl;
//    std::cout << "Newton's method" << std::endl;
//    std::cout << "x = " << ans1 << std::endl;
//    std::cout << "Iterations: " << iterNewton << std::endl;
//    std::cout << "Simple iterations method" << std::endl;
//    std::cout << "x = " << ans2 << std::endl;
//    std::cout << "Iterations: " << iterSimple << std::endl;
//    return 0;
//}
