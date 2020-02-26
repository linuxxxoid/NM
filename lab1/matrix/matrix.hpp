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


class Matrix {
public:
    Matrix();
    Matrix(int32_t n, int32_t m);
    
    std::vector<double>& operator[](const int32_t index);
    const std::vector<double>& operator[](const int index) const;
    int32_t getRow() const;
    std::vector<std::vector<double>> getMatrix();
    void getKoefficients(std::vector<std::vector<double>>&);
    
    bool isSquareMatrix() const;
    bool isTridiagonalMatrix() const;
    bool isSimmetricalMatrix() const;
    
    ~Matrix();
private:
    std::vector<std::vector<double>> _matrix;
    int32_t row, column;
};

#endif /* matrix_hpp */
