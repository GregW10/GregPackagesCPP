#include "simsup.hpp"

using ld_pair = std::pair<long double, long double>;

#define NUM_SIMS 33

gtd::bod_0f jupiter{189813*BILLION*BILLION*10*10*10*10, 69'911'000, {}, {}};

unsigned int factor = 8;
long double dt = 1.0l/factor;
uint64_t iterations = 100*factor;
uint64_t num_frames = 2500;

const long double starting_angle = 0.0l; // rad
// const long double end_bd = 3'000.0l; // kg/m^3
const long double angle_step = gtd::PI/32; // rad

const long double bulk_density = 500;

long double comet_rad = 750.0l;
long double body_rad = 70.0l;

gtd::vector3D<long double> comet_pos{-414139744.3484770l, 277323640.2236369l, -1231468367.968793l};
gtd::vector3D<long double> comet_vel{3491.263406628809l, -6314.208154956334l, 11582.30048080498l};

const long double spacing = 0.1l;
long double body_mass{1};

const long double cor = 1.0l;

gtd::vec3 orientation{0, 0, 1};
const gtd::vec3 omega;

std::pair<uint64_t, long double> bods_rad() {
    auto [_sys, crad] = gtd::sys::hcp_comet(comet_rad, comet_pos, comet_vel, spacing, body_mass, body_rad, cor,
                                            orientation, omega, true, true);
    return {_sys.num_bodies(), crad};
}

bool btrk_f(const gtd::bod_0f &_b) {
    return _b.mass() < BILLION*BILLION;
}

void output_data(std::ostream& os,
                 long double real_rad,
                 const std::vector<std::tuple<long double, gtd::vec3, std::string>>& tups,
                 const ld_pair *seps) {
    os << "time_elapsed,target_comet_radius,actual_comet_radius,bulk_density,body_mass,body_spacing,body_cor,w_x,w_y,w_z\n" <<
    dt*iterations*num_frames << ',' << comet_rad << ',' << real_rad << ',' << bulk_density << ',' << body_mass << ',' << spacing << ',' << cor << ',' <<
    omega.get_x() << ',' << omega.get_y() << ',' << omega.get_z() << '\n' <<
    "angle,or_x,or_y,or_z,mean_sep,mean_com_dist,,,,\n";
    for (const auto &[angle, ornt, _] : tups)
        os << angle << ',' << ornt.get_x() << ',' << ornt.get_y() << ',' << ornt.get_z() << ',' << seps->first << ',' <<
        seps++->second << ",,,," << '\n';
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr,"Inavlid number of command-line arguments.\nUsage: ./bdsim <output_csv_path> [opt:num_sims].\n");
        return 1;
    }
    int64_t num_sims = NUM_SIMS;
    if (argc == 3)
        num_sims = gtd::to_int<int64_t>(*(argv + 2));
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (!max_threads)
        max_threads = num_sims;
    std::cout << "Number of parallel threads that will be run: " << max_threads << std::endl;
    std::pair<uint64_t, long double> pair = bods_rad();
    body_mass = gtd::comet_body_mass(pair.second, bulk_density, pair.first);
    printf("Number of bodies in comet: %" PRIu64"\nComet effective radius: %.30Lf\nBody mass: %.30Lf\n",
           pair.first, pair.second, body_mass);
    std::vector<std::thread> threads;
    threads.reserve(max_threads);
    long double angle = starting_angle;
    auto btrk_func = [](const gtd::bod_0f &_b){return _b < BILLION*BILLION;};
    uint64_t counter = 0;
    std::vector<std::tuple<long double, gtd::vec3, std::string>> vals;//, long double, long double>> vals;
    vals.reserve(num_sims);
    const gtd::vec3 perp = gtd::vec_ops::cross(comet_pos, comet_vel);
    while (counter < num_sims) {
        vals.emplace_back(angle, comet_pos.copy().rodrigues_rotate(perp, angle),
                          std::to_string(angle).insert(0, "Angle"));//, 0, 0);
        angle = starting_angle + angle_step*(++counter);
    }
    ld_pair *seps = new ld_pair[num_sims]{};
    ld_pair *sptr = seps;
    counter = 1;
    for (auto &[theta, ornt, path] : vals) {
        threads.emplace_back(&gtd::run_sim<long double, decltype(&btrk_f),
        decltype(&gtd::system<long double, long double, long double, false, false, 3, 0, 0, false>::hcp_comet<false>),
        long double, gtd::vec3, gtd::vec3, long double, long double, long double, long double, gtd::vec3, gtd::vec3,
        bool, bool, long double, uint64_t, const char *>,
                             sptr++, path.c_str(), "log", 30, std::ref(jupiter), dt, iterations, num_frames, &btrk_f,
                   &gtd::system<long double, long double, long double, false, false, 3, 0, 0, false>::hcp_comet<false>,
                             comet_rad, comet_pos, comet_vel, spacing, body_mass, body_rad, cor,
                             ornt, omega, true, true, 1.0l, 1l, "M1:1,D1:1,T1:1");
        std::cout << "----------------\nSimulation " << counter << '/' << num_sims << "\nAngle: " << theta <<
        " rad\nOrientation: " << ornt << " \n----------------" << std::endl;
        if (!(counter++ % max_threads)) {
            for (std::thread &_t : threads)
                _t.join();
            threads.clear();
        }
    }
    if (!threads.empty())
        for (std::thread &_t : threads)
            _t.join();
    if (std::ofstream csv{*(argv + 1), std::ios_base::out | std::ios_base::trunc}; csv) {
        csv.precision(21);
        output_data(csv, pair.second, vals, seps);
        csv.close();
    } else {
        std::cout << "Error opening output file. Writing data to standard output instead:\n----------------\n";
        std::cout.precision(21);
        output_data(std::cout, pair.second, vals, seps);
        std::cout.flush();
    }
    if (FILE *fp = fopen("sep_pair_data.psep", "wb"); fp != nullptr) {
        fwrite(seps, sizeof(ld_pair), num_sims, fp);
        fclose(fp);
    }
    delete [] seps;
    return 0;
}
