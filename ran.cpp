#include "gregsys.hpp"
#include "gregastro.hpp"

const gtd::vec3 pos{0, 0, 0};
const gtd::vec3 vel{0, 0, 0};
const uint64_t num = 100;
const long double bounding_r = 1000.0l;
const long double b_mass = 26500000000;
const long double b_rad = 75.0l;
const long double restc_f = 1.0l;

uint64_t num_reps = 15000;

int main() {
	gtd::system<long double, long double, long double, true, false, 3, 0, 0, false> sys =
            gtd::system<long double, long double, long double, true, false, 3, 0, 0, false>::
                    random_comet(pos, vel, num, bounding_r, b_mass, b_rad, restc_f, gtd::sys::leapfrog_kdk, 0.25l,
                                 0, 0.01);
    time_t id = time(nullptr);
    sys.set_iterations(1000);
    sys.set_timestep(0.1l);
    for (auto& bod : sys) {
        bod.set_restitution(0.5l);
    }
    gtd::image_dimensions dims{1000, 1000};
    gtd::asc_0f asc{dims.x, dims.y};
    asc.add_system(sys);
    asc.set_num_decor_stars(0);
    gtd::cam cam{dims};
    cam.set_position({0, -3000, 0});
    cam.set_direction({0, 1, 0});
    asc.follow_camera(&cam);
    gtd::star_t star{1, 1, {0, 0, 100'000}, {}, 1, 1};
    asc.add_body(star);
    gtd::star_t star2{1, 1, {0, 100'000}, {}, 1, 1};
    asc.add_body(star2);
    gtd::star_t star3{1, 1, {100'000, 0}, {}, 1, 1};
    asc.add_body(star3);
    gtd::star_t star4{1, 1, {0, 0, -100'000}, {}, 1, 1};
    asc.add_body(star4);
    gtd::star_t star5{1, 1, {0, -100'000}, {}, 1, 1};
    asc.add_body(star5);
    gtd::star_t star6{1, 1, {-100'000, 0}, {}, 1, 1};
    asc.add_body(star6);
    uint64_t counter = 0;
    long long num_zeros;
    long double base = log10l(num_reps);
    gtd::String path;
    gtd::String npath;
    uint64_t num_reps_half = num_reps/6;
    sys.to_nsys(gtd::String{".nsys"}.append_front(id).c_str());
    unsigned int i;
    for (i = 1; i <= num_reps_half; ++i) {
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
    }
    for (auto& bod : sys) {
        bod.set_restitution(1.0l);
    }
    for (; i <= num_reps; ++i) {
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
    }
    asc.render();
    asc.write("ImageZ.bmp"); // guaranteed to always be last
    sys.to_nsys("Data_FileZ.nsys");
	return 0;
}
