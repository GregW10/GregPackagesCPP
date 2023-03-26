#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>

int main(int argc, char **argv) {
    if (argc < 4 || argc > 5) {
        std::cerr << "Invalid number of command-line arguments.\n"
                     "Expected usage: ./mass <comet_radius> <bulk_density> <num_bods> [opt:precision]\n";
        return 1;
    }
    long double radius;
    long double bulk_density;
    unsigned long long num_bods;
    int precision = 10;
    try {
        radius = std::stold(*(argv + 1));
        bulk_density = std::stold(*(argv + 2));
        num_bods = std::stoull(*(argv + 3));
        if (argc == 5)
            precision = std::stoi(*(argv + 4));
    } catch (const std::invalid_argument& e) {
        std::cout << "Error! Conversion could not be performed.\nwhat(): " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Mass of each body: " << std::setprecision(precision) <<
    (((4.0l/3)*M_PI*radius*radius*radius)*bulk_density)/num_bods << " kg" << std::endl;
    return 0;
}
