//
//  matrix.cpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#include "matrix.hpp"

Matrix::Matrix() {
    row = 0;
    column = 0;
}

Matrix::Matrix(int32_t n, int32_t m) {
    try {
        if (n < 0 || m < 0) {
            throw 1;
        }
        row = n;
        column = m;
        _matrix.assign(row, std::vector<double>(column, 0));
    }
    catch (int32_t i) {
        std::cout << "Error №" << i << ": size of matrix must be > 0!" << std::endl;
    }
}

std::ostream& operator<<(std::ostream& os, const Matrix& m) {
    os << "Matrix " << m.row << "x" << m.column << ": " << std::endl;
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.column; ++j) {
            os << '\t';
            os << m[i][j];
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

std::vector<double>& Matrix::operator[](const int32_t index) {
    return _matrix[index];
}

const std::vector<double>& Matrix::operator[](const int index) const {
    return _matrix[index];
}

const Matrix operator*(const Matrix& left, double rightNumber){
    Matrix result = left;
    for(int i = 0; i < result.getRow(); ++i){
        for(int j = 0; j < result.getColumn(); ++j){
            result[i][j] *= rightNumber;
        }
    }
    return result;
}

const Matrix operator*(double leftNumber, const Matrix& right) {
    return right * leftNumber;
}

//const Matrix operator*(double leftNumber, Matrix& right) {
//    return right * leftNumber;
//}

const Matrix operator*(const Matrix& left, const Matrix& right) {
    if (left.getColumn() != right.getRow()) {
        throw "Wrong sizes of matrixes to multiply!";
    }
    Matrix result(left.getRow(), right.getColumn());
    for (int i = 0; i < result.getRow(); ++i) {
        for (int j = 0; j < result.getColumn(); ++j) {
            for (int m = 0; m < left.getColumn(); ++m) {
                result[i][j] += left[i][m] * right[m][j];
            }
        }
    }
    return result;
}

const std::vector<double> operator*(const Matrix& left, const std::vector<double>& right) {
    if (left.getColumn() != right.size()) {
        throw "Wrong sizes of matrixes to multiply!";
    }
    std::vector<double> res(left.getRow(), 0.0);
    for(int i = 0; i < left.getRow(); ++i){
        for(int j = 0; j < left.getColumn(); ++j){
           res[i] += left[i][j] * right[j];
        }
    }
    return res;
}

const Matrix operator-(const Matrix& left, const Matrix& right){
    if(left.getRow() != right.getRow() || left.getColumn() != right.getColumn()){
        throw "Wrong minus! Sizes of matrix not equal!";
    }
    Matrix ans(left.getRow(), left.getColumn());
    for(int i = 0; i < ans.getRow(); ++i){
        for(int j = 0; j < ans.getColumn(); ++j){
            ans[i][j] = left._matrix[i][j] - right._matrix[i][j];
        }
    }
    return ans;
}

int32_t Matrix::getRow() const {
    return row;
}

int32_t Matrix::getColumn() const {
    return column;
}

std::vector<std::vector<double>> Matrix::getMatrix() {
    return _matrix;
}

void Matrix::getKoefficients(std::vector<std::vector<double>>& koef) {
    try {
        if (!isSquareMatrix()) {
            throw 2;
        }
        koef.assign(row, std::vector<double>(3, 0.0));

        for (int i = 0; i < row; ++i) {
            koef[i][0] = i - 1 < 0 ? 0.0 : _matrix[i][i - 1];
            koef[i][1] = _matrix[i][i];
            koef[i][2] = i + 1 == column ? 0.0 : _matrix[i][i + 1];
        }
    } catch (int32_t i) {
        std::cout << "Error №" << i << ": it's not a square matrix!" << std::endl;
        return;
    }
}


void Matrix::unitary() {
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < column; ++j) {
            _matrix[i][j] = i != j ? 0.0 : 1.0;
        }
    }
}

void Matrix::transpose() {
    std::vector<std::vector<double>> tmp(column, std::vector<double>(row));
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < column; ++j) {
            tmp[j][i] = _matrix[i][j];
        }
    }
    _matrix = tmp;
}


bool Matrix::isSquareMatrix() const {
    return row == column;
}

bool Matrix::isTridiagonalMatrix() const {
    try {
        if(!isSquareMatrix()) {
            throw 2;
        }
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < column; ++j) {
                if ((abs(i - j) > 1) && _matrix[i][j]){
                    return false;
                }
            }
        }
        return true;
    }
    catch (int32_t i) {
        std::cout << "Error №" << i << ": it's not a square matrix!" << std::endl;
        return false;
    }
}

bool Matrix::isSimmetricalMatrix() const{
    try {
        if(!isSquareMatrix()) {
            throw 2;
        }
        for (int i = 0; i < row; ++i){
            for (int j = i + 1; j < column; ++j) {
                if(_matrix[i][j] != _matrix[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }
    catch (int32_t i) {
        std::cout << "Error №" << i << ": it's not a square matrix!" << std::endl;
        return false;
    }
}


Matrix::~Matrix() {}
