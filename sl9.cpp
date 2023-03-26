// #define GREGSYS_MERGERS

#include "simsup.hpp"

gtd::bod_0f jupiter{189813*BILLION*BILLION*10*10*10*10, 69'911'000, {}, {}};
// gtd::bod_0f b2{500'000'000.0l, 10, {2970, 2222, 1111}, {}};
// gtd::bod_0f b3{500'000'000.0l, 10, {3000, 2222, 1141}, {}};
// gtd::bod_0f b4{500'000'000.0l, 10, {3000, 2222, 1101}, {}};

long double dt = 0.125;
uint64_t iterations = 800;
size_t num_reps = 2521;

long double comet_rad = 750.0l;
long double b_rad = 75.0l;
long double b_mass = 3646105641.631895542l;
long double b_sep = 0.1l;

gtd::vector3D<long double> pos{-414139744.3484770l, 277323640.2236369l, -1231468367.968793l};
gtd::vector3D<long double> vel{3491.263406628809l, -6314.208154956334l, 11582.30048080498l};

int main(int argc, char **argv) {
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    time_t id = time(nullptr);
    gtd::String starting_time_str{gtd::get_date_and_time()};
    /*                              comet rad  comet pos comet vel.    sep   b.m.  b.rad.  r_coeff */
    gtd::system<long double, long double, long double, true, false, 3, 0, 0, false> sys =
            gtd::system<long double, long double, long double, true, false, 3, 0, 0, false>::
                    hcp_comet(comet_rad, pos, vel, b_sep, b_mass, b_rad, 1, true);
    sys.add_body(jupiter);
    sys.set_iterations(iterations);
    sys.set_timestep(dt);
#ifdef GREGSYS_MERGERS
    sys.set_min_tot_com_mom(BILLION*BILLION*BILLION*BILLION*BILLION*BILLION);
#endif
    gtd::String nsys_path;
    nsys_path.append_back(id).append_back(".nsys");
    sys.to_nsys(nsys_path.c_str());
    gtd::image_dimensions dims = {500, 500};
    gtd::asc_0f asc{dims.x, dims.y};
    gtd::star_t star{1, 1, {0, 0, 2'000'000'000}, {}, 1, 1};
    //asc.add_star(star);
    gtd::star_t star2{1, 1, {-2'000'000'000, 0, 0}, {}, 1, 1};
    //asc.add_star(star2);
    gtd::star_t star3{1, 1, {0, -2'000'000'000, 0}, {}, 1, 1};
    //asc.add_star(star3);
    gtd::star_t star4{1, 1, {2'000'000'000, 0}, {}, 1, 1};
    //asc.add_star(star4);
    gtd::star_t star5{1, 1, {0, 2'000'000'000, 0}, {}, 1, 1};
    //asc.add_star(star5);
    gtd::star_t star6{1, 1, {0, 0, -2'000'000'000}, {}, 1, 1};
    //asc.add_star(star6);
    gtd::cam cam;
    cam.set_image_dimensions(dims);
    // cam.set_position({3000, 1900, 1111});
    cam.set_direction({0, 1, 0});
    asc.follow_camera(&cam);
    asc.set_num_decor_stars(0);
    gtd::btrk_0f btrk{sys, [](const gtd::bod_0f& b){return b.mass() < pow(10, 20);}};
    long long num_zeros;
    long double base = log10l(num_reps);
    asc.add_system(sys);
    asc.set_body_clr(sys.back().get_id(), gtd::colors::aqua); // back() is Jupiter
    std::cout << "Number of bodies: " << sys.num_bodies() << std::endl;
    gtd::String path;
    std::cout << "Com vel: " << btrk.com_vel() << std::endl;
    gtd::vec3 dir = gtd::vec_ops::cross(btrk.com(), btrk.com_vel());
    gtd::String npath;
    gtd::vec3 com;
    gtd::vec3 comv;
    long double time_per_step = iterations*dt;
    std::cout << "Time: 0 seconds" << std::endl;
    // long double dist = 15'000;
    long double distf = 6;
    // gtd::vec3 offset = ;
    gtd::vec3 dir_hat = dir.unit_vector();
    gtd::vec3 changing_dir;
    cam.set_direction(dir);
    gtd::star_t the_one{1, 1, gtd::vec3{1.69424*powl(10,7), -6.44198*powl(10,7), 6.63999*powl(10,7)}.unit_vector()*3'000'000'000, {}, 1, 1};
    if (argc == 2) {
        asc.add_star(the_one);
    } else {
        asc.add_star(star);
        asc.add_star(star2);
        asc.add_star(star3);
        asc.add_star(star4);
        asc.add_star(star5);
        asc.add_star(star6);
    }
    for (unsigned int i = 1; i <= num_reps; ++i) {
        // asc.add_system(sys);
        com = btrk.com();
        comv = btrk.com_vel();
        if (argc == 2) {
            dir_hat = (jupiter.pos() - com).unit_vector();
            cam.set_direction(dir_hat);
            cam.set_position({com - dir_hat*distf*btrk.avg_dist_to(com)});
        } else if (argc == 3) {
            dir_hat = gtd::vec_ops::cross(btrk.com(), btrk.com_vel()).unit_vector();
            dir_hat.rodrigues_rotate(comv, -gtd::PI/2);
            cam.set_direction(dir_hat);
            cam.set_position({com - dir_hat*distf*btrk.avg_dist_to(com)});
        } else {
            cam.set_position({com - dir_hat*distf*btrk.avg_dist_to(com)});
        }
        // cam.set_direction({btrk.com() - cam.position()});
        asc.render();
        path = "Image";
        npath = "Data_File";
        num_zeros = base - floorl(log10l(i));
        for (unsigned int j = 0; j < num_zeros; ++j) {
            path.append_back("0");
            npath.append_back("0");
        }
        path.append_back(i).append_back(".bmp");
        npath.append_back(i).append_back(".nsys");
        asc.write(path.c_str());
        sys.to_nsys(npath.c_str());
        sys.evolve(gtd::sys::leapfrog_kdk);
        std::cout << "Image " << +i << '/' << num_reps << " written." << std::endl;
        path.clear();
        // asc.clear_bodies();
        std::cout << "COM: " << com << "\nCOM vel.: " << comv << std::endl;
        std::cout << "COM distance: " << (com - jupiter.pos()).magnitude() <<
        "\nCOM vel. mag.: " << comv.magnitude() << std::endl;
        std::cout << "Time: " << time_per_step*i << " seconds" << std::endl;
    }
    asc.render();
    asc.write("ImageZ.bmp"); // guaranteed to always be last
    sys.to_nsys("Data_FileZ.nsys");
    gtd::String ending_time_str{gtd::get_date_and_time()};
    // std::cout << sys << std::endl;
    std::chrono::time_point<std::chrono::high_resolution_clock> finish = std::chrono::high_resolution_clock::now();
#ifdef __linux__
    std::cout << "Completed rendering and simulation at: " << gtd::get_date_and_time() << std::endl;
    std::cout << "Total time taken to complete: " << (finish - start).count()/BILLION << " seconds." << std::endl;
    path.clear();
    path = "Comet_Video_";
    path.append_back(gtd::get_date_and_time()).append_back("_nsys").append_back(id);
    gtd::String path_mp4 = path;
    path.append_back(".y4m");
    path_mp4.append_back(".mp4");
    std::cout << "\nGenerating .y4m video...\n" << std::endl;
    std::system(("ym -ptsd -o " + path).c_str());
    std::cout << "\nGenerating .mp4 video...\n" << std::endl;
    std::system(("ffmpeg -y -i "+ path + " -c:v libx264 -crf 18 -preset veryslow -pix_fmt yuv420p -movflags +faststart "
                                        + path_mp4).c_str());
    std::cout << "\nMoving .mp4 to bucket with gsutil...\n" << std::endl;
    std::system(path_mp4.append_front("gsutil mv ").append_back(" gs://nbod_bucket/").c_str());
    std::cout << "\nCopying .nsys to bucket with gsutil...\n" << std::endl;
    std::system(nsys_path.append_front("gsutil cp ").append_back(" gs://nbod_bucket/").c_str());
    gtd::update_log("nbod_bucket", "logs.txt", 60, id, "nbodn", starting_time_str.c_str(), ending_time_str.c_str(),
                    sys.num_bodies(), dt, iterations, num_reps + 1,
                    sys.num_bodies() - 1, comet_rad, b_rad, b_mass, b_sep, pos, vel);
#endif
    return 0;
}
