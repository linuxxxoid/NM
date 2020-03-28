//
//  main.cpp
//  lab3
//
//  Created by Лина Вельтман on 24/03/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

//Interpolation Polynomial

double F(double x) {
    return sqrt(x);
}

double LagrangeMethod(std::vector<double>& x, double point_x) {
    double res = 0;
    for (int i = 0; i < x.size(); ++i) {
        double f = F(x[i]);
        for (int j = 0; j < x.size(); ++j) {
            if (i != j) {
                f *= (point_x - x[j]) / (x[i] - x[j]);
            }
        }
        res += f;
    }
    return abs(F(point_x) - res);
}

double w(std::vector<double>& x, double k, double ii) {
    double res = 1;
    for (int i = 0; i < ii + 1; ++i) {
        if (i != k) {
            res *= (x[k] - x[i]);
        }
    }
    return res;
}

std::vector<double> SeparatedDifferences(std::vector<double>& x) {
    std::vector<double> f(x.size()), A(x.size());
    for (int i = 0; i < x.size(); ++i) {
        f[i] = F(x[i]);
    }
    
    for (int i = 0; i < x.size(); ++i) {
        for (int k = 0; k <= i; ++k) {
            A[i] += f[k] / w(x, k, i);
        }
    }
    return A;
}

double NewtonMethod(std::vector<double>& x, double point_x) {
    int size = (int) x.size();
    std::vector<double> A = SeparatedDifferences(x);
    double res = 0, difference = 1;
    for (int i = 0; i < size; ++i) {
        res += A[i] * difference;
        difference *= point_x - x[i];
    }
    return abs(F(point_x) - res);
}

void InterpolationPolynomial() {
    double x = 3.0;
    std::vector<double> a({ 0, 1.7, 3.4, 5.1 });
    std::vector<double> b({ 0, 1.7, 4.0, 5.1 });

    std::cout << "Lagrange Method:" << std::endl;
    std::cout << "a)\t" <<  LagrangeMethod(a, x) << std::endl;
    std::cout << "b)\t" <<  LagrangeMethod(b, x) << std::endl;

    std::cout << "Newton Method:" << std::endl;
    std::cout << "a)\t" <<  NewtonMethod(a, x) << std::endl;
    std::cout << "b)\t" <<  NewtonMethod(b, x) << std::endl;
}

//Cube Spline
int FindInterval(std::vector<double>& x, double point_x);

double FuncS(double a, double b, double c, double d, double x) {
    return a + b * x + c * x * x + d * x * x * x;
}

std::vector<double> getA(std::vector<double>& f) {
    int n = (int) f.size();
    std::vector<double> a(n, 0);
    for (int i = 0; i < n; ++i) {
        a[i] = f[i];
    }
    return a;
}

std::vector<double> getB(std::vector<double>& f, std::vector<double>& h, std::vector<double>& c) {
    int n = (int) h.size();
    std::vector<double> b(n, 0);
    for (int i = 0; i < n - 1; ++i) {
        b[i] = (f[i + 1] - f[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
    }
    b[n - 1] = (f[n] - f[n - 1]) / h[n - 1] - 2 * h[n - 1] * c[n - 1] / 3;
    return b;
}

std::vector<double> getC(std::vector<double>& points, std::vector<double>& values, std::vector<double>& h, std::vector<double>& a) {
    int size = (int) points.size();
    std::vector<double> sub_a(size - 2);
    std::vector<double> sub_b(size - 2);
    std::vector<double> sub_c(size - 2);
    std::vector<double> sub_d(size - 2);
    for (int i = 0; i < size - 3; i++) {
        sub_b[i] = h[i];
        sub_a[i + 1] = h[i];
    }
    for (int i = 0; i < size - 2; i++) {
        sub_c[i] = 2 * (h[i + 1] + h[i]);
    }
    for (int i = 1; i < size - 1; i++) {
        sub_d[i - 1] = 3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]);
    }
    
    int n = size - 2;
    std::vector<double> res(n + 1);
    res[0] = 0;
    for (int i = 1; i < n; i++) {
        sub_c[i] -= sub_a[i] * sub_b[i - 1] / sub_c[i - 1];
        sub_d[i] -= sub_a[i] * sub_d[i - 1] / sub_c[i - 1];

    }
    res[n] = sub_d[n - 1] / sub_c[n - 1];
    for (int i = n - 1; i > 0; --i) {
        res[i] = (sub_d[i - 1] - sub_b[i - 1] * res[i + 1]) / sub_c[i - 1];
    }
    return res;
}

std::vector<double> getD(std::vector<double>& h, std::vector<double>& c) {
    int n = (int) h.size();
    std::vector<double> d(n, 0);
    for (int i = 0; i < n - 1; ++i) {
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
    d[n - 1] = -c[c.size() - 1] / (3 * h[n - 1]);
    return d;
}

void WriteToFile(std::vector<double>& points, std::vector<double>& values, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
    std::ofstream file;
    file.open("log2.txt");
    if (file.is_open()) {
        for (int i = 0; i < points.size(); ++i) {
            file << points[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < values.size(); ++i) {
            file << values[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < a.size(); ++i) {
            file << a[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < b.size(); ++i) {
            file << b[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < c.size(); ++i) {
            file << c[i] << " ";
        }
        file << std::endl;
        for (int i = 0; i < d.size(); ++i) {
            file << d[i] << " ";
        }
        file << std::endl;
        file.close();
    }
}


void CubeSplineInterpolation() {
//    std::vector<double> points({0.0, 1.0, 2.0, 3.0, 4.0});
//    std::vector<double> values({0.0, 1.8415, 2.9093, 3.1411, 3.2432});
    //    double point_x = 1.5;

    std::vector<double> points({0.0, 1.7, 3.4, 5.1, 6.8});
    std::vector<double> values({0.0, 1.3038, 1.8439, 2.2583, 2.6077});
    double point_x = 3.0;

    std::vector<double> h(points.size() - 1);
    for (int i = 0; i < h.size(); ++i) {
        h[i] = points[i + 1] - points[i];
    }
    std::vector<double> a = getA(values);

    std::vector<double> c = getC(points, values, h, a);
    std::vector<double> b = getB(values, h, c);
    std::vector<double> d = getD(h, c);
    
    
    int i = FindInterval(points, point_x);
    double res = FuncS(a[i], b[i], c[i], d[i], point_x - points[i]);
    std::cout << "Cube Spline Interpolation" << std::endl;
    std::cout << "Value: " << res << std::endl;
    
    WriteToFile(points, values, a, b, c, d);
    
}

//

double SumOfSquaredErrors(std::vector<double>& f, std::vector<double>& y) {
    double res = 0.0;
    for (int i = 0; i < f.size(); ++i) {
        res += (f[i] - y[i]) * (f[i] - y[i]);
    }
    return res;
}

//Differencials

int FindInterval(std::vector<double>& x, double point_x) {
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] <= point_x && point_x <= x[i + 1]) {
            return i;
        }
    }
    return -1;
}

double DiffLeft(std::vector<double>& x, std::vector<double>& y, double point_x) {
    double i = FindInterval(x, point_x);
    if (i == -1) {
        std::cout << "Error!\n";
        return -1;
    }
    return (y[i + 1] - y [i]) / (x[i + 1] - x[i]);
}

double DiffRight(std::vector<double>& x, std::vector<double>& y, double point_x) {
    int i = -1;
    for (i = (int) x.size() - 2; i >= 0; --i) {
        if (x[i + 1] >= point_x && x[i] <= point_x) {
            break;
        }
    }
    if (i == -1) {
        std::cout << "Error!\n";
        return -1;
    }
    return (y[i + 1] - y [i]) / (x[i + 1] - x[i]);
}

double FirstDiff(std::vector<double>& x, std::vector<double>& y, double point_x) {
    double i = FindInterval(x, point_x);
    if (i == -1) {
        std::cout << "Error!\n";
        return -1;
    }
    double first = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    double second = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    double third = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    return ((first - second) / (x[i + 2] - x[i])) * (2 * point_x - x[i] - x[i + 1]) + third;
}

double SecondDiff(std::vector<double>& x, std::vector<double>& y, double point_x) {
    double i = FindInterval(x, point_x);
    if (i == -1) {
        std::cout << "Error!\n";
        return -1;
    }
    double first = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
    double second = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    return 2 * (first - second) / (x[i + 2] - x[i]);
}

void DiffMethod() {
    double point_x = 0.2;
    std::vector<double> x({ -0.2, 0.0, 0.2, 0.4, 0.6 });
    std::vector<double> y({ 1.7722, 1.5708, 1.3694, 1.1593, 0.9273 });

    std::cout << "First derivative with first order accuracy" << std::endl;
    std::cout << "Left-side derivative: " << DiffLeft(x, y, point_x) << std::endl;
    std::cout << "Right-side derivative: " << DiffRight(x, y, point_x) << std::endl;
    std::cout << "First derivative with second order accuracy: " << FirstDiff(x, y, point_x) << std::endl;
    std::cout << "Second derivative: " << SecondDiff(x, y, point_x) << std::endl;
}

//Integration

double FuncToIntegrate(double x) {
    return 1 / (3 * x * x + 4 * x + 2);
//    return x / ((3 * x + 4) * (3 * x + 4));
}

double RectangleMethod(std::vector<double>& x, double h) {
    double sum = 0;
    for (int i = 0; i < x.size() - 1; ++i) {
        sum += FuncToIntegrate((x[i] + x[i + 1]) / 2);
    }
    return h * sum;
}

double TrapezeMethod(std::vector<double>& y, double h) {
    double sum = 0;
    int n = (int) y.size() - 1;
    for (int i = 1; i < n; ++i) {
        sum += y[i];
    }
    return h * (y[0] / 2 + sum + y[n] / 2);
}

double SimpsonMethod(std::vector<double>& y, double h) {
    //четный, нечетный
    double even = 0, odd = 0;
    int n = (int) y.size() - 1;
    for (int i = 1; i < n; ++i) {
        if (i & 1) {
            odd += y[i];
        } else {
            even += y[i];
        }
    }
    return (h / 3) * (y[0] + 4 * odd + 2 * even + y[n]);
}

void RungeRombergRichardsonMethod(std::vector<std::vector<double>>& res, double h1, double h2) {
    double anal_val = 1.8574186872187473;
    double k = 2;
    std::pair<double, double> rec = {abs(res[0][0] - res[1][0]) / (k * k - 1), abs(res[0][0] - anal_val) / (k * k - 1)};
    std::pair<double, double> trap = {abs(res[0][1] - res[1][1]) / (k * k - 1), abs(res[0][1] - anal_val) / (k * k - 1)};
    std::pair<double, double> simp = {abs(res[0][2] - res[1][2]) / (k * k - 1), abs(res[0][2] - anal_val) / (k * k - 1)};
    std::cout << "Runge-Romberg-Richardson Method" << std::endl;
    std::cout << "Rectangle error: "<< rec.first << "\t" << rec.second << std::endl;
    std::cout << "Trapeze error: " << trap.first << "\t" << trap.second << std::endl;
    std::cout << "Simpson error: " << simp.first << "\t" << simp.second << std::endl;
}

void get_points(std::vector<double>& points, double x0, double xk, double step) {
    points.clear();
    double val = x0;
    while (val < xk + step) {
        points.emplace_back(val);
        val += step;
    }
}

void get_values(std::vector<double>& values, std::vector<double>& points) {
    values.clear();
    for (int i = 0; i < points.size(); ++i) {
        values.emplace_back(FuncToIntegrate(points[i]));
    }
}

void IntegrationMethod() {
    std::vector<std::vector<double>> results(2);
    double x0 = -2, xk = 2, h1 = 1.0, h2 = 0.5;
//    double x0 = -1, xk = 1, h1 = 0.5, h2 = 0.25;
    
    std::vector<double> points, values;
    get_points(points, x0, xk, h1);
    get_values(values, points);

    double rec, trap, simp;
    
    rec = RectangleMethod(points, h1);
    results[0].emplace_back(rec);
    std::cout << "Rectangle method with step = " << h1 << std::endl;
    std::cout << rec << std::endl;
    
    trap = TrapezeMethod(values, h1);
    results[0].emplace_back(trap);
    std::cout << "Trapeze method with step = " << h1 << std::endl;
    std::cout << trap << std::endl;
    
    simp = SimpsonMethod(values, h1);
    results[0].emplace_back(simp);
    std::cout << "Simpson method with step = " << h1 << std::endl;
    std::cout << simp << std::endl;
    
    
    get_points(points, x0, xk, h2);
    get_values(values, points);
    
    rec = RectangleMethod(points, h2);
    results[1].emplace_back(rec);
    std::cout << "Rectangle method with step = " << h2 << std::endl;
    std::cout << rec << std::endl;
    
    trap = TrapezeMethod(values, h2);
    results[1].emplace_back(trap);
    std::cout << "Trapeze method with step = " << h2 << std::endl;
    std::cout << trap << std::endl;
    
    simp = SimpsonMethod(values, h2);
    results[1].emplace_back(simp);
    std::cout << "Simpson method with step = " << h2 << std::endl;
    std::cout << simp << std::endl;
    
    RungeRombergRichardsonMethod(results, h1, h2);
}


void menu() {
    std::cout << "=================" << std::endl;
    std::cout << "Choose:" << std::endl;
    std::cout << "1 - Construction of interpolation polynomials." << std::endl;
    std::cout << "2 - Building a cube spline." << std::endl;
    std::cout << "3 - Finding approximate polynomials." << std::endl;
    std::cout << "4 - Computations of derivatives." << std::endl;
    std::cout << "5 - Calculation of a definite integral." << std::endl;
    std::cout << "0 - Exit" << std::endl;
    std::cout << "=================" << std::endl;
}

int main(int argc, const char * argv[]) {
    while (1) {
        menu();
        int choice;
        std::cin >> choice;
        if (choice == 1) {
            InterpolationPolynomial();
        }
        else if (choice == 2) {
            CubeSplineInterpolation();
        }
        else if (choice == 3) {
            
        }
        else if (choice == 4) {
            DiffMethod();
        }
        else if (choice == 5) {
            IntegrationMethod();
        }
        else if (choice == 0) {
            std::cout << "End!\n";
            break;
        }
        else {
            std::cout << "Input error! Try again\n";
        }
    }
    return 0;
}
