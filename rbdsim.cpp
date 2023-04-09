#define RANDOM_COMET

#include "simsup.hpp"

using ld_pair = std::pair<long double, long double>;
using sys_t = gtd::system<long double, long double, long double, false, false, 3, 0, 0, false>;

#define NUM_SIMS 30

gtd::bod_0f jupiter{189813*BILLION*BILLION*10*10*10*10, 69'911'000, {}, {}};

unsigned int factor = 8;
long double dt = 1.0l/factor;
uint64_t iterations = 100*factor;
uint64_t num_frames = 2500;

const long double starting_bd = 100.0l; // kg/m^3
// const long double end_bd = 3'000.0l; // kg/m^3
const long double bd_step = 100.0l; // kg/m^3

/* * * * * * * * * * * * * * * * * * * */
const uint64_t num_bods = 907;
const long double bounding_rad = 1'600.0l;
const long double min_cor = 0.5l;
const long double max_cor = 0.95l;
const long double n_exp = 3.0l;
const long double d_scale = 1200.0l;
const long double ev_dt = 0.0625l/8.0l;
const uint64_t probing_iters = 100'000;
const long double dist_tol = 0.0l;

const long double avg_pf = 0.6568027311608282l;
/* * * * * * * * * * * * * * * * * * * */

// long double comet_rad = 750.0l;
long double body_rad = 70.0l;

auto b_massf = [](const long double &rho){return (SPHERE_VOLUME(body_rad)*rho)/avg_pf;};

gtd::vector3D<long double> comet_pos{-414139744.3484770l, 277323640.2236369l, -1231468367.968793l};
gtd::vector3D<long double> comet_vel{3491.263406628809l, -6314.208154956334l, 11582.30048080498l};

const long double body_temp_mass = 1.0l;

const long double cor = 1.0l;

const gtd::vec3 omega;

// std::pair<uint64_t, long double> bods_rad() {
//     auto [_sys, crad] = gtd::sys::hcp_comet(comet_rad, comet_pos, comet_vel, spacing, body_temp_mass, body_rad, cor,
//                                             orientation, omega, true, true);
//     return {_sys.num_bodies(), crad};
// }

bool btrk_f(const gtd::bod_0f &_b) {
    return _b.mass() < BILLION*BILLION;
}

void output_data(std::ostream& os,
                 gtd::com_props<long double> *_p,
                 uint64_t num_sims) {
    os << "time_elapsed,target_pf,body_cor,w_x,w_y,w_z\n" <<
    dt*iterations*num_frames << ',' << avg_pf << ',' << cor << ',' <<
    omega.get_x() << ',' << omega.get_y() << ',' << omega.get_z() << '\n' <<
    "bulk_density,effective_radius,body_mass,pf,mean_sep,mean_com_dist\n";
    uint64_t counter = 0;
    while (counter++ < num_sims)
        os << _p->_bd << ',' << _p->_effr << ',' << _p->_b_mass << ',' << _p->_pf << ',' << _p->_mu_sep << ',' <<
        _p++->_mu_csep << '\n';
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
    // std::pair<uint64_t, long double> pair = bods_rad();
    printf("Number of bodies in comet: %" PRIu64"\n", num_bods);
    std::vector<std::thread> threads;
    threads.reserve(max_threads);
    long double bd = starting_bd;
    auto btrk_func = [](const gtd::bod_0f &_b){return _b < BILLION*BILLION;};
    //                     bulk density body mass   directory     avg. sep.   sep. from com.
    // std::vector<com_props<long double>> vals;//, long double, long double>> vals;
    // vals.reserve(num_sims);
    // while (counter++ < num_sims) {
    //     vals.emplace_back(bd, 0, 0, 0, 0, gtd::comet_body_mass(pair.second, bd, pair.first),
    //                       std::to_string(bd).insert(0, "BulkDensity"));//, 0, 0);
    //     bd += bd_step;
    // }
    gtd::com_props<long double> *props = new gtd::com_props<long double>[num_sims]{};
    gtd::com_props<long double> *pptr = props;
    // ld_pair *seps = new ld_pair[num_sims]{};
    // ld_pair *sptr = seps;
    std::vector<std::string> paths;
    paths.reserve(max_threads);
    uint64_t counter = 1;
    while (counter <= num_sims) { // yes, creating lambdas willy-nilly is bad, but this is a throwaway program
        auto lam = [bd](sys_t &_sys, const long double &_pf){return gtd::adjust_bd(_sys, body_rad, _pf, bd);};
        threads.emplace_back(&gtd::run_sim<long double, decltype(lam), decltype(&btrk_f),
        std::tuple<sys_t, long double, long double> (*)(const gtd::vec3&, const gtd::vec3&, const gtd::vec3&, uint64_t,
                const long double&, const long double&, const long double&, long double, long double, long double, int,
                long double, long double, long double, uint64_t, long double, uint64_t, long double, const char*),
        gtd::vec3, gtd::vec3, gtd::vec3, uint64_t, long double, long double, long double, long double, long double, long double,
        int, long double, long double, long double, uint64_t, long double, uint64_t, long double, const char *>,
                             pptr++, paths.emplace_back(std::to_string(bd)).insert(0, "BulkDensity").c_str(), "log", 30,
                             std::ref(jupiter), dt, iterations, num_frames,
                             lam,
                             &btrk_f, (std::tuple<sys_t, long double, long double> (*)(const gtd::vec3&, const gtd::vec3&, const gtd::vec3&, uint64_t,
                                                                                       const long double&, const long double&, const long double&, long double, long double, long double, int,
                                                                                       long double, long double, long double, uint64_t, long double, uint64_t, long double, const char*)) &sys_t::random_comet<false>,
                             /*****************************************************************************************/
                             comet_pos, comet_vel, omega, num_bods, bounding_rad, b_massf(bd), body_rad, cor, n_exp, d_scale,
                             sys_t::leapfrog_kdk, min_cor, max_cor, ev_dt, probing_iters, dist_tol, iterations, dt,
                             "M1:1,D1:1,T1:1");
        std::cout << "----------------\nSimulation " << counter << '/' << num_sims << "\nBulk density: " << bd <<
        " kg/m^3\nBody starting mass: " << b_massf(bd) << " kg\n----------------" << std::endl;
        // pptr++->_bd = bd;
        bd += bd_step;
        if (!(counter++ % max_threads)) {
            for (std::thread &_t : threads)
                _t.join();
            threads.clear();
            paths.clear();
        }
    }
    if (!threads.empty())
        for (std::thread &_t : threads)
            _t.join();
    if (std::ofstream csv{*(argv + 1), std::ios_base::out | std::ios_base::trunc}; csv) {
        csv.precision(21);
        output_data(csv, props, num_sims);
        csv.close();
    } else {
        std::cout << "Error opening output file. Writing data to standard output instead:\n----------------\n";
        std::cout.precision(21);
        output_data(std::cout, props, num_sims);
        std::cout.flush();
    }
    if (FILE *fp = fopen("sep_data.props", "wb"); fp != nullptr) {
        fwrite(props, sizeof(gtd::com_props<long double>), num_sims, fp);
        fclose(fp);
    }
    delete [] props;
    return 0;
}
