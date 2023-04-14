#include "simsup.hpp"

bool funky = 1;

int main(int argc, char **argv) {
    if (argc != 6) {
        std::cerr << "Invalid number of command-line arguments."
                     "\nUsage: ./renderall <nsys_prefix> <nsys_freq> <width> <height> <cam_option>\n";
        return 1;
    }
    const char *prefix = *(argv + 1);
    unsigned long long freq = std::stoull(*(argv + 2)); // will throw if cannot convert
    unsigned long width = std::stoul(*(argv + 3));
    unsigned long height = std::stoul(*(argv + 4));
    unsigned char cam_option = std::stoull(*(argv + 5));
    long int size;
    if ((size = pathconf(".", _PC_PATH_MAX)) == -1)
        size = 4095;
    char *curr_dir = new char[size + 1]{};
    if (getcwd(curr_dir, size) == nullptr) {
        delete [] curr_dir;
        std::cerr << "Error obtaining current working directory.\n";
        abort();
    }
    std::vector<std::string> paths;
    namespace fs = std::filesystem;
    for (const fs::directory_entry &entry : fs::directory_iterator(std::string{curr_dir}))
        if (entry.path().string().ends_with(".nsys") && entry.path().filename().string().starts_with(prefix))
            paths.emplace_back(entry.path().filename().string());
    std::sort(paths.begin(), paths.end());
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    gtd::String starting_time_str{gtd::get_date_and_time()};
    gtd::system<long double, long double, long double, false, false, 3, 0, 0, false> sys{paths.front().c_str()};
    gtd::image_dimensions dims = {(unsigned int) width, (unsigned int) height};
    gtd::asc_0f asc{dims.x, dims.y};
    gtd::star_t star{1, 1, {0, 0, 2'000'000'000}, {}, 1, 1};
    gtd::star_t star2{1, 1, {-2'000'000'000, 0, 0}, {}, 1, 1};
    gtd::star_t star3{1, 1, {0, -2'000'000'000, 0}, {}, 1, 1};
    gtd::star_t star4{1, 1, {2'000'000'000, 0}, {}, 1, 1};
    gtd::star_t star5{1, 1, {0, 2'000'000'000, 0}, {}, 1, 1};
    gtd::star_t star6{1, 1, {0, 0, -2'000'000'000}, {}, 1, 1};
    gtd::cam cam;
    cam.set_image_dimensions(dims);
    asc.follow_camera(&cam);
    asc.set_num_decor_stars(0);
    gtd::btrk_0f btrk{sys, [](const gtd::bod_0f& b){return b.mass() < pow(10, 20);}};
    asc.set_body_col_gen(gtd::vibrant_col_gen(128, 70));
    asc.add_system(sys);
    asc.set_body_clr(sys.back().get_id(), gtd::colors::aqua); // back() is Jupiter
    if (funky)
        for (const auto &[bid, clr] : gtd::col_dist(btrk, gtd::colors::red, gtd::colors::green, gtd::colors::blue)) {
            asc.set_body_clr(bid, clr);
        }
    std::cout << "Number of bodies: " << sys.num_bodies() << std::endl;
    std::cout << "Comet number of bodies: " << btrk.num_bods() << std::endl;
    gtd::String path;
    gtd::vec3 dir = gtd::vec_ops::cross(btrk.com(), btrk.com_vel());
    gtd::vec3 com;
    gtd::vec3 comv;
    long double distf = 6;
    gtd::vec3 dir_hat = dir.unit_vector();
    cam.set_direction(dir);
    gtd::star_t the_one{1, 1, gtd::vec3{1.69424*powl(10,7), -6.44198*powl(10,7), 6.63999*powl(10,7)}.unit_vector()*3'000'000'000, {}, 1, 1};
    if (cam_option == 2) {
        asc.add_star(the_one);
    } else {
        asc.add_star(star);
        asc.add_star(star2);
        asc.add_star(star3);
        asc.add_star(star4);
        asc.add_star(star5);
        asc.add_star(star6);
    }
    gtd::bod_0f jupiter = sys.back();
    unsigned long long counter = 0;
    std::vector<std::string>::size_type first_digit_index = gtd::strlen_c(prefix);
    // for (const char &_c : paths.front()) {
    //     if (_c >= 48 || _c <= 57)
    //         break;
    //     ++first_digit_index;
    // }
    gtd::bod_0f *ptr;
    for (const std::string &_path : paths) {
        if (counter % freq) {
            ++counter;
            continue;
        }
        gtd::system<long double, long double, long double, false, false, 3, 0, 0, false> cpy_from{_path.c_str()};
        ptr = const_cast<gtd::bod_0f*>(&sys.front());
        for (const gtd::bod_0f &_b : cpy_from)
            ptr++->update(_b.pos(), _b.vel());
        com = btrk.com();
        comv = btrk.com_vel();
        if (cam_option == 2) {
            dir_hat = (jupiter.pos() - com).unit_vector();
            cam.set_direction(dir_hat);
            cam.set_position({com - dir_hat*distf*btrk.mean_dist_to(com)});
        } else if (cam_option == 3) {
            dir_hat = gtd::vec_ops::cross(btrk.com(), btrk.com_vel()).unit_vector();
            dir_hat.rodrigues_rotate(comv, -gtd::PI/2);
            cam.set_direction(dir_hat);
            cam.set_position({com - dir_hat*distf*btrk.mean_dist_to(com)});
        } else {
            cam.set_position({com - dir_hat*distf*btrk.mean_dist_to(com)});
        }
        asc.render();
        path = "Image";
        for (const char &c : _path.substr(first_digit_index)) {
            if (c == '.')
                break;
            path.push_back(c);
        }
        path.append_back(".bmp");
        asc.write(path.c_str());
        std::cout << "Image " << ++counter << " written\r";
        std::cout.flush();
        path.clear();
    }
    gtd::String ending_time_str{gtd::get_date_and_time()};
    std::chrono::time_point<std::chrono::high_resolution_clock> finish = std::chrono::high_resolution_clock::now();
    std::cout << "\nTime taken: " <<
    std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start).count()/BILLION << " seconds" << std::endl;
    return 0;
}
