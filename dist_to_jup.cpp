#include "simsup.hpp"

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Invalid number of command-line arguments.\n";
        return 1;
    }
    gtd::system<long double, long double, long double, true, false, 3, 0, 0, false> sys{*(argv + 1)};
    const gtd::bod_0f& jupiter = sys.back();
    gtd::btrk_0f btrk{sys, [](const gtd::bod_0f& _b){return _b.mass() < BILLION*BILLION;}};
    std::cout.precision(30);
    const gtd::vec3 com = btrk.com();
    std::cout << "Distance: " << (com - jupiter.pos()).magnitude() << std::endl;
    std::cout << "Mean separation: " << btrk.mean_sep() << std::endl;
    std::cout << "Mean distance to COM: " << btrk.mean_dist_to(com) << std::endl;
    return 0;
}