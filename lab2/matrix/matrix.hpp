//
//  matrix.hpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#ifndef matrix_hpp
#define matrix_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>


class Matrix {
public:
    Matrix();
    Matrix(int32_t n, int32_t m);
    
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
    std::vector<double>& operator[](const int32_t index);
    const std::vector<double>& operator[](const int index) const;
    friend const Matrix operator*(const Matrix& left, double rightNumber);
    friend const Matrix operator*(double leftNumber, const Matrix& right);
    friend const std::vector<double> operator*(const Matrix& left, const std::vector<double>& right);
    //friend const Matrix operator*(double leftNumber, Matrix& right);
    friend const Matrix operator*(const Matrix& left, const Matrix& right);
    friend const Matrix operator-(const Matrix& left, const Matrix& right);

    int32_t getRow() const;
    int32_t getColumn() const;
    std::vector<std::vector<double>> getMatrix();
    void getKoefficients(std::vector<std::vector<double>>&);
    
    void unitary();
    void transpose();
    
    bool isSquareMatrix() const;
    bool isTridiagonalMatrix() const;
    bool isSimmetricalMatrix() const;
    
    ~Matrix();
private:
    std::vector<std::vector<double>> _matrix;
    int32_t row, column;
};

#endif /* matrix_hpp */
