#include <iostream>
#include <string>
#include <cmath>

#define PI 3.14159265358979323846264338327950288419716939937510582097494459230

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Invalid number of command-line arguments."
                     "\nUsage: ./rad <desired_comet_radius> <number_of_bodies> <packing_fraction>\n";
        return 1;
    }
    auto sphv = [](auto radius){return (4.0l/3.0l)*PI*radius*radius*radius;}; // I like lambdas ok??
    auto rad = [](auto volume){return cbrtl((3*volume)/(4*PI));};
    long double c_rad = std::stold(*(argv + 1));
    unsigned long long num = std::stoull(*(argv + 2));
    long double pf = std::stold(*(argv + 3));
    printf("\033[33mRadius of each body: \033[34m%.30Lf\n", rad(sphv(c_rad)*pf/num));
    return 0;
}
