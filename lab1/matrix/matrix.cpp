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

std::vector<double>& Matrix::operator[](const int32_t index) {
    return _matrix[index];
}

const std::vector<double>& Matrix::operator[](const int index) const {
    return _matrix[index];
}

int32_t Matrix::getRow() const {
    return row;
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
