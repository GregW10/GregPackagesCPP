//
// Created by Gregor Hartl Watters on 20/03/2023.
//

#include "gregsys.hpp"
#include <iomanip>

template <gtd::isNumWrapper T>
void print_reduced_sys(const gtd::system<T, T, T, true, false, 3, 0, 0, false>& sys) {
    uint64_t counter = 1;
    for (const auto& _b : sys) {
        std::cout << "Body " << counter++ << ":\nMass = " << _b.mass() << "\nRadius = " << _b.rad() << "\nPosition = "
        << _b.pos() << "\nVelocity = " << _b.vel() << "\n--------------------\n";
    }
    std::cout.flush();
}

template <gtd::isNumWrapper T, uint64_t rF>
void print_info(const gtd::body_tracker<T, T, T, rF>& btrk) {
    gtd::vector3D<T> com = btrk.com();
    gtd::vector3D<T> comv = btrk.com_vel();
    uint64_t num_bods = btrk.num_bods();
    long double mean_sep = 0;
    const gtd::body<T, T, T, rF> *const *outer = btrk.front();
    const gtd::body<T, T, T, rF> *const *inner;
    uint64_t oc = 0;
    uint64_t ic;
    while (oc < num_bods) {
        inner = outer + 1;
        ic = ++oc;
        while (ic++ < num_bods) {
            mean_sep += gtd::vec_ops::distance((*inner)->pos(), (*outer)->pos());
            ++inner;
        }
        ++outer;
    }
    mean_sep /= (num_bods*(num_bods - 1))/2;
    std::cout << "--------------------\nFiltered bodies info:\n--------------------\n";
    std::cout << "Number of bodies: " << num_bods << "\nTotal mass: " << btrk.mass() << '\n' <<
    "COM pos.: " << com << "\nCOM vel.: " << comv << '\n' <<
    "COM pos. mag.: " << com.magnitude() << "\nCOM vel. mag.: " << comv.magnitude() << "\nAverage distance to COM: " <<
    btrk.mean_dist_to(com) << "\nMean sep. between bodies: " << mean_sep << std::endl;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "Error: invalid number of command-line arguments.\nUsage: scat [nsys_filename] "
                     "[opt:mass_threshold] [opt:precision]\n";
        return 1;
    }
    if (argc == 4)
        std::cout << std::setprecision(std::stoi(*(argv + 3)));
    std::cout << "System bodies:\n-------------\n\n";
    try {
        gtd::system<long double, long double, long double, true, false, 3, 0, 0, false> sys{*(argv + 1)};
        print_reduced_sys(sys);
        if (argc >= 3) {
            long double mass_cutoff = std::stold(*(argv + 2)); // throws std::invalid_argument in case of failure
            gtd::btrk_0f btrk{sys, [&mass_cutoff](const gtd::bod_0f& bod){return bod.mass() < mass_cutoff;}};
            print_info(btrk);
        }
    } catch (const std::exception& e) {
        try {
            gtd::system<double, double, double, true, false, 3, 0, 0, false> sys{*(argv + 1)};
            print_reduced_sys(sys);
            if (argc >= 3) {
                double mass_cutoff = std::stod(*(argv + 2)); // throws std::invalid_argument in case of failure
                gtd::body_tracker<double, double, double, 0>
                btrk{sys, [&mass_cutoff](const gtd::body<double, double, double, 0>& bod){return bod.mass() < mass_cutoff;}};
                print_info(btrk);
            }
        } catch (const std::exception& e) {
            try {
                gtd::system<float, float, float, true, false, 3, 0, 0, false> sys{*(argv + 1)};
                print_reduced_sys(sys);
                if (argc >= 3) {
                    float mass_cutoff = std::stod(*(argv + 2)); // throws std::invalid_argument in case of failure
                    gtd::body_tracker<float, float, float, 0>
                    btrk{sys, [&mass_cutoff](const gtd::body<float, float, float, 0>& bod)
                    {return bod.mass() < mass_cutoff;}};
                    print_info(btrk);
                }
            } catch (const std::exception& e) {
                std::cerr << "Error loading .nsys file, what(): " << e.what() << '\n';
                return 1;
            }
        }
    }
    return 0;
}
