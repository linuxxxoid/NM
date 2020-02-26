//
//  main.cpp
//  lab1
//
//  Created by Лина Вельтман on 26/02/2020.
//  Copyright © 2020 Lina Veltman. All rights reserved.
//

#include "lab1_1.hpp"
#include "lab1_2.hpp"
#include "lab1_3.hpp"
#include "lab1_4.hpp"
#include "lab1_5.hpp"


void help() {
    std::cout << "Enter parameter: " << std::endl;
    std::cout << 1 << " - to solve SLAU by LU separating matrix with Gauss method." << std::endl;
    std::cout << "Find determinant, back matrix, L and U matrix. " << std::endl;
    std::cout << 2 << " - to solve SLAU by Thomas algorithm (progonka). " << std::endl;
    std::cout << 3 << " - to solve SLAU by simple ittearation and Zeidel methods." << std::endl;
    std::cout << "Analyze the number of iterations which are required to achieve a given accuracy. " << std::endl;
    std::cout << 4 << "Enter parameter: " << std::endl;
    std::cout << 5 << "Enter parameter: " << std::endl;

}

int main(int argc, const char * argv[]) {

    std::cout << "LABORATORY WORK №1" << std::endl;
    std::cout << "Variant №7" << std::endl;
    help();
    
    int choice;
    while (std::cin >> choice) {
        if (choice == 1) {
//            if (!lab1_1()) {
//                std::cout << "Something wrong! Try again!\n";
//            }
        }
        else if (choice == 2) {
            if (!lab1_2()) {
                std::cout << "Something wrong! Try again!\n";
            }
        }
        else if (choice == 3) {
//            if (!lab1_3()) {
//                std::cout << "Something wrong! Try again!\n";
//            }
        }
        else if (choice == 4) {
//            if (!lab1_4()) {
//                std::cout << "Something wrong! Try again!\n";
//            }
        }
        else if (choice == 5) {
//            if (!lab1_5()) {
//                std::cout << "Something wrong! Try again!\n";
//            }
        }
        else {
            std::cout << "Error input!\n";
            help();
        }
    }
    return 0;
}
