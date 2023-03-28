#include "gregsys.hpp"

template <gtd::isNumWrapper T>
void output_energies(const gtd::system<T, T, T, true, false, 3, 0, 0, false>& sys) {
    auto num = sys.num_bodies();
    long double ke = 0;
    long double pe = 0;
    long double *bpe = new long double[num]{};
    long double *cpy = bpe;
    const gtd::body<T, T, T, 0> *bod_o = &sys.front();
    const gtd::body<T, T, T, 0> *bod_i;
    uint64_t i = 0;
    uint64_t j;
    long double gr;
    long double G = sys.G_val();
    while (i < num) {
        ke += 0.5l*bod_o->mass()*(bod_o->vel()*bod_o->vel());
        for (j = i + 1, bod_i = bod_o + 1; j < num; ++j, ++bod_i) {
            gr = -G/(bod_i->pos() - bod_o->pos()).magnitude();
            *bpe += bod_i->mass()*gr;
            *(bpe + j - i) += bod_o->mass()*gr;
        }
        pe += *bpe++;
        ++i;
        ++bod_o;
    }
    delete [] cpy;
    std::cout << "--------------------\nKinetic Energy: " << ke << "\nPotential Energy: " << pe << "\nTotal Energy: " <<
    (ke + pe) << "\n--------------------\n";
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Error: Invalid number of command line arguments.\n";
        return 1;
    }
    for (char **ptr = argv + 1; --argc > 0; ++ptr) {
        std::cout << "--------------------\n" << *ptr << ":\n";
        try {
            gtd::system<long double, long double, long double, true, false, 3, 0, 0, false> sys{*ptr};
            output_energies(sys);
            continue;
        } catch (const gtd::nbody_error& err) {
            std::cerr << "long double NOT matched.\nwhat(): " << err.what();
        }
        try {
            gtd::system<double, double, double, true, false, 3, 0, 0, false> sys{*ptr};
            output_energies(sys);
            continue;
        } catch (const gtd::nbody_error& err) {
            std::cerr << "double NOT matched.\nwhat(): " << err.what();
        }
        try {
            gtd::system<float, float, float, true, false, 3, 0, 0, false> sys{*ptr};
            output_energies(sys);
            continue;
        } catch (const gtd::nbody_error& err) {
            std::cerr << "float NOT matched.\n";
            std::cerr << "Error loading nsys file \"" << *ptr << "\". " << "what(): " << err.what() << '\n';
            return 1;
        }
        ++ptr;
    }
    return 0;
}
