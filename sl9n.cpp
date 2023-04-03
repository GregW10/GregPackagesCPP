// #define GREGSYS_MERGERS

#include "simsup.hpp"
#include <iomanip>

gtd::bod_0f jupiter{189813*BILLION*BILLION*10*10*10*10, 69'911'000, {}, {}};
// gtd::bod_0f b2{500'000'000.0l, 10, {2970, 2222, 1111}, {}};
// gtd::bod_0f b3{500'000'000.0l, 10, {3000, 2222, 1141}, {}};
// gtd::bod_0f b4{500'000'000.0l, 10, {3000, 2222, 1101}, {}};

unsigned int factor = 8;
long double dt = 1.0l/factor;
uint64_t iterations = 100*factor;
size_t num_reps = 2521;

long double comet_rad = 750.0l;
long double b_rad = 75.0l;
long double b_mass = 132513706.19408596907l;
long double b_sep = 0.1l;

gtd::vector3D<long double> pos{-414139744.3484770l, 277323640.2236369l, -1231468367.968793l};
gtd::vector3D<long double> vel{3491.263406628809l, -6314.208154956334l, 11582.30048080498l};

int main(int argc, char **argv) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    time_t id = time(nullptr);
    gtd::String starting_time_str{gtd::get_date_and_time()};
    /*                              comet rad  comet pos comet vel.    sep   b.m.  b.rad.  r_coeff */
    auto [sys, crad] =
            gtd::system<long double, long double, long double, true, false, 3, 0, 0, false>::
                    hcp_comet<false>(comet_rad, pos, vel, b_sep, b_mass, b_rad, 1, true);
    printf("Comet effective radius: %.30Lf\n", crad);
    sys.add_body(jupiter, false);
    sys.set_iterations(iterations);
    sys.set_timestep(dt);
#ifdef GREGSYS_MERGERS
    sys.set_min_tot_com_mom(BILLION*BILLION*BILLION*BILLION*BILLION*BILLION);
#endif
    gtd::String nsys_path;
    nsys_path.append_back(id).append_back(".nsys");
    sys.to_nsys(nsys_path.c_str());
    gtd::btrk_0f btrk{sys, [](const gtd::bod_0f& b){return b.mass() < pow(10, 20);}};
    long long num_zeros;
    long double base = log10l(num_reps);
    std::cout << "Number of bodies: " << sys.num_bodies() << std::endl;
    std::cout << "Comet number of bodies: " << btrk.num_bods() << std::endl;
    std::cout << "Com vel: " << btrk.com_vel() << std::endl;
    gtd::String npath;
    gtd::vec3 com;
    gtd::vec3 comv;
    long double time_per_step = iterations*dt;
    std::cout << "Time: 0 seconds" << std::endl;
    long double distf = 4;
    gtd::vec3 dir_hat1 = gtd::vec_ops::cross(btrk.com(), btrk.com_vel()).unit_vector();
    gtd::vec3 dir_hat2;
    gtd::vec3 dir_hat3;
    gtd::vec3 cam_pos1;
    gtd::vec3 cam_pos2;
    gtd::vec3 cam_pos3;
    std::ofstream cpar;
    gtd::String cpar_path;
    int prec = 20;
    if (argc == 2) {
        prec = std::stoi(*(argv + 1));
    }
    for (unsigned int i = 1; i <= num_reps; ++i) {
        com = btrk.com();
        comv = btrk.com_vel();
        npath = "Data_File";
        cpar_path = "Camera1_Params";
        num_zeros = base - floorl(log10l(i));
        for (unsigned int j = 0; j < num_zeros; ++j) {
            npath.append_back("0");
            cpar_path.append_back("0");
        }
        npath.append_back(i).append_back(".nsys");
        cpar_path.append_back(i).append_back(".cpar");
        cam_pos1 = com - dir_hat1*distf*btrk.avg_dist_to(com);
        dir_hat2 = (jupiter.pos() - com).unit_vector();
        cam_pos2 = com - dir_hat2*distf*btrk.avg_dist_to(com);
        dir_hat3 = gtd::vec_ops::cross(btrk.com(), btrk.com_vel()).unit_vector();
        dir_hat3.rodrigues_rotate(comv, -gtd::PI/2);
        cam_pos3 = com - dir_hat3*distf*btrk.avg_dist_to(com);
        cpar.open(cpar_path.c_str(), std::ios_base::out | std::ios_base::trunc);
        cpar << std::setprecision(prec) <<
        "Image dimensions = 4000 4000\n\n" << "Camera:\nPosition = " << cam_pos1.get_x() << ' ' <<
        cam_pos1.get_y() << ' ' << cam_pos1.get_z() << "\nDirection = " << dir_hat1.get_x() << ' ' <<
        dir_hat1.get_y() << ' ' << dir_hat1.get_z() << "\nRotation (rad) = 0\nImage distance = 2\n\nStar:\nMass = 1\n"
                                                     "Radius = 1\nPosition = 0 0 2000000000\nVelocity = 0 0 0\n"
                                                     "Luminosity = 1\nPoint sources = 1\n\n"
                                                     "Star:\nMass = 1\n"
                                                     "Radius = 1\nPosition = 0 2000000000 0\nVelocity = 0 0 0\n"
                                                     "Luminosity = 1\nPoint sources = 1\n\n"
                                                     "Star:\nMass = 1\n"
                                                     "Radius = 1\nPosition = 2000000000 0 0\nVelocity = 0 0 0\n"
                                                     "Luminosity = 1\nPoint sources = 1\n\n"
                                                     "Star:\nMass = 1\n"
                                                     "Radius = 1\nPosition = 0 0 -2000000000\nVelocity = 0 0 0\n"
                                                     "Luminosity = 1\nPoint sources = 1\n\n"
                                                     "Star:\nMass = 1\n"
                                                     "Radius = 1\nPosition = 0 -2000000000 0\nVelocity = 0 0 0\n"
                                                     "Luminosity = 1\nPoint sources = 1\n\n"
                                                     "Star:\nMass = 1\n"
                                                     "Radius = 1\nPosition = -2000000000 0 0\nVelocity = 0 0 0\n"
                                                     "Luminosity = 1\nPoint sources = 1\n\n";
        cpar.close();
        cpar_path[6] = '2';
        cpar.open(cpar_path.c_str(), std::ios_base::out | std::ios_base::trunc);
        cpar << "Image dimensions = 4000 4000\n\n" << "Camera:\nPosition = " << cam_pos2.get_x() << ' ' <<
        cam_pos2.get_y() << ' ' << cam_pos2.get_z() << "\nDirection = " << dir_hat2.get_x() << ' ' <<
        dir_hat2.get_y() << ' ' << dir_hat2.get_z() << "\nRotation (rad) = 0\nImage distance = 2\n\nStar:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 0 2000000000\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 2000000000 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 2000000000 0 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 0 -2000000000\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 -2000000000 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = -2000000000 0 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n";
        cpar.close();
        cpar_path[6] = '3';
        cpar.open(cpar_path.c_str(), std::ios_base::out | std::ios_base::trunc);
        cpar << "Image dimensions = 4000 4000\n\n" << "Camera:\nPosition = " << cam_pos3.get_x() << ' ' <<
        cam_pos3.get_y() << ' ' << cam_pos3.get_z() << "\nDirection = " << dir_hat3.get_x() << ' ' <<
        dir_hat3.get_y() << ' ' << dir_hat3.get_z() << "\nRotation (rad) = 0\nImage distance = 2\n\nStar:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 0 2000000000\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 2000000000 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 2000000000 0 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 0 -2000000000\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = 0 -2000000000 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n"
                                                            "Star:\nMass = 1\n"
                                                            "Radius = 1\nPosition = -2000000000 0 0\nVelocity = 0 0 0\n"
                                                            "Luminosity = 1\nPoint sources = 1\n\n";
        cpar.close();
        sys.to_nsys(npath.c_str());
        sys.evolve(gtd::sys::leapfrog_kdk);
        // std::cout << sys.num_bodies() << std::endl;
        // return 1;
        std::cout << "Iteration " << +i << '/' << num_reps << " completed." << std::endl;
        // asc.clear_bodies();
        std::cout << "COM: " << com << "\nCOM vel.: " << comv << std::endl;
        std::cout << "COM distance: " << (com - jupiter.pos()).magnitude() <<
        "\nCOM vel. mag.: " << comv.magnitude() << std::endl;
        std::cout << "Time: " << time_per_step*i << " seconds" << std::endl;
    }
    sys.to_nsys("Data_FileZ.nsys");
    gtd::String ending_time_str{gtd::get_date_and_time()};
    // std::cout << sys << std::endl;
    std::chrono::time_point<std::chrono::high_resolution_clock> finish = std::chrono::high_resolution_clock::now();
#ifdef __linux__
    std::cout << "Completed rendering and simulation at: " << gtd::get_date_and_time() << std::endl;
    std::cout << "Total time taken to complete: " << (finish - start).count()/BILLION << " seconds." << std::endl;
    std::cout << "\nCopying .nsys to bucket with gsutil...\n" << std::endl;
    std::system(nsys_path.append_front("gsutil cp ").append_back(" gs://nbod_bucket/").c_str());
    gtd::update_log("nbod_bucket", "logs.txt", 60, id, "nbod", starting_time_str.c_str(), ending_time_str.c_str(),
                    sys.num_bodies(), dt, iterations, num_reps + 1,
                    btrk.num_bods(), crad, b_rad, b_mass, b_sep, pos, vel);
#endif
    return 0;
}
