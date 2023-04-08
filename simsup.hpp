#ifndef GREG_SIMSUP_H // simulation support
#define GREG_SIMSUP_H

#ifdef _WIN32
#error "This header file is for UNIX systems only.\n"
#endif

#include "gregsys.hpp"
#include "gregastro.hpp"
#include <sstream>

namespace gtd {
    gtd::String generation_id(const char *bucket, const char *fname) {
        if (!bucket || !fname || !*bucket || !*fname) {
            return false;
        }
        gtd::String gen;
        gtd::String cmd = "gsutil stat gs://";
        cmd.append_back(bucket).push_back('/').append_back(fname);
        FILE *fp;
        if ((fp = popen(cmd.c_str(), "r")) == nullptr) {
            std::cerr << "Error with popen().\n";
            abort();
        }
        int c;
        while ((c = fgetc(fp)) != EOF) {
            if (c == '\n') {
                if (gen.contains("Generation")) {
                    gen.erase_chars(0, gen.find(':') + 1);
                    gen.strip(" \t\n");
                    break;
                }
                gen.erase_chars();
                continue;
            }
            gen.push_back((char) c);
        }
        if (int retval = pclose(fp)) {
            if (retval == -1) {
                std::cerr << "Error with popen().\n";
                abort();
            }
            std::cerr << "Error associated with command execution.\n";
            abort();
        }
        return gen;
    }
    bool update_text(const char *bucket, const char *fname, const char *txt_to_append, unsigned int wait_time = 10) {
        if (!bucket || !fname || !txt_to_append || !*bucket || !*fname || !*txt_to_append || !wait_time) {
            return false;
        }
        gtd::String gen;
        gtd::String download_cmd = "gsutil cp gs://";
        download_cmd.append_back(bucket).push_back('/').append_back(fname).append_back(" .");
        gtd::String upload_cmd = "gsutil -h x-goog-if-generation-match:";
        std::ofstream f;
        while (true) {
            gen = generation_id(bucket, fname);
            if (std::system(download_cmd.c_str())) {
                std::cerr << "Error download original file with gsutil.\n";
                abort();
            }
            f.open(fname, std::ios_base::app);
            f << txt_to_append;
            f.close();
            if (generation_id(bucket, fname) == gen)
                break;
            sleep(wait_time);
        }
        upload_cmd.append_back(gen).
        append_back(" mv ").append_back(fname).append_back(" gs://").append_back(bucket).push_back('/');
        if (std::system(upload_cmd.c_str())) {
            std::cerr << "Error uploading file to bucket.\n";
            abort();
        }
        return true;
    }
    void update_log(const char *bucket, const char *log_fname, unsigned int wait_time, time_t id, const char *server,
                    const char *starting_time, const char *ending_time,
                    uint64_t num_bods, long double dt, uint64_t steps_per_frame, uint64_t frames, uint64_t c_num_bods,
                    long double pf, long double c_rad, long double c_prad, long double c_pmass, long double c_psep,
                    const vec3& c_pos, const vec3& c_vel, bool random_comet = false,
                    const vec3& orientation = vec3::zero, const vec3& omega = vec3::zero, const char *mass_units = "kg",
                    const char *dist_units = "m", const char *time_units = "s") {
        if (!bucket || !log_fname || !*bucket || !*log_fname || !wait_time)
            return;
        long int size;
        if ((size = pathconf(".", _PC_PATH_MAX)) == -1)
            size = 4095;
        char *dir = new char[size + 1]{};
        if (getcwd(dir, size) == nullptr) {
            delete [] dir;
            std::cerr << "Error obtaining current working directory.\n";
            abort();
        }
        std::ostringstream oss; // to store entire text to be written to file
        oss << "Simulation ID: " << id << "\nServer: " << server << "\nDirectory: " << dir;
        delete [] dir;
        oss.precision(30);
        // if (!random_comet) {
        oss << "\nStart time: " << starting_time << "\nEnd time: " << ending_time << "\nNumber of bodies: " <<
        num_bods << "\nTime-step: " << dt << ' ' << time_units <<
        "\nSteps per frame: " << steps_per_frame << "\nFrames: " << frames << "\nRandom comet: " <<
        std::boolalpha << random_comet << "\nPacking fraction: " << pf << "\nComet number of bodies: " <<
        c_num_bods << "\nComet radius: " << c_rad << ' ' << dist_units << "\nComet particle radius: " << c_prad <<
        ' ' << dist_units << "\nComet particle mass: " << c_pmass << ' ' << mass_units <<
        "\nComet particle separation: " << c_psep << ' ' << dist_units << "\nComet position: " << c_pos << ' ' <<
        dist_units << "\nComet velocity: " << c_vel << ' ' << dist_units << '/' << time_units <<
        "\nComet orientation: " << orientation << "\nComet angular velocity: " << omega << "rad/s" <<
        "\n----------------------------------------\n";
        // } else {
        //     oss << "\nStart time: " << starting_time << "\nEnd time: " << ending_time << "\nNumber of bodies: " <<
        //         num_bods << "\nTime-step: " << dt << ' ' << time_units <<
        //         "\nSteps per frame: " << steps_per_frame << "\nFrames: " << frames << "\nComet number of bodies: " <<
        //         c_num_bods << "\nComet radius: N/A\nComet particle radius: " << c_prad <<
        //         ' ' << dist_units << "\nComet particle mass: " << c_pmass << ' ' << mass_units <<
        //         "\nComet particle separation: N/A\nComet position: " << c_pos << ' ' <<
        //         dist_units << "\nComet velocity: " << c_vel << ' ' << dist_units << '/' << time_units <<
        //         "\n----------------------------------------\n";
        // }
        update_text(bucket, log_fname, oss.str().c_str(), wait_time);
    }
    template <isNumWrapper T, isNumWrapper M>
    T comet_bd(const T& comet_radius, const M& comet_mass) {
        if (comet_radius < 0 || comet_mass < 0)
            throw std::invalid_argument{"Error: invalid parameter/s.\n"};
        return comet_mass/SPHERE_VOLUME(comet_radius);
    }
    template <isNumWrapper T>
    T comet_body_mass(const T& comet_radius, const T& bulk_density, uint64_t num_bods) {
        if (comet_radius < 0 || bulk_density < 0)
            throw std::invalid_argument{"Error: invalid parameter/s.\n"};
        return (SPHERE_VOLUME(comet_radius)*bulk_density)/num_bods;
    }
    template <typename sys_t>
    long double adjust_bd(sys_t &sys, long double b_rad, long double pf, long double new_bd) { // adjust bulk density of comet
        if (!sys.num_bodies())
            return std::numeric_limits<long double>::quiet_NaN();
        long double b_mass = (SPHERE_VOLUME(b_rad)*new_bd)/pf;
        for (auto &_b : sys)
            _b.set_mass(b_mass);
        return b_mass;
    }
    auto vibrant_col_gen(long double _min_avg, long double _min_sd) {
        return [_min_avg, _min_sd](){
            static std::mt19937_64 engine{(uint64_t) time(nullptr)};
            std::uniform_int_distribution<unsigned char> dist{0, 255};
            gtd::color retc;
            do {
                retc.b = dist(engine);
                retc.g = dist(engine);
                retc.r = dist(engine);
            } while (sd(retc.b, retc.g, retc.r) < _min_sd || mean_avg(retc.b, retc.g, retc.r) < _min_avg);
            return retc;
        };
    }
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t mF>
    std::map<uint64_t, color> col_dist(const body_tracker<M, R, T, mF> &btrk,
                                       const color &_h, // hot
                                       const color &_i, // intermediate
                                       const color &_c) /* cold */ {
        if (!btrk.num_bods())
	        return {};
	    vector3D<T> _com = btrk.com();
        T _furthest_distance{};
        T _distance;
        std::vector<std::pair<uint64_t, T>> distances;
	    for (const body<M, R, T, mF> *const &_btr : btrk) {
            if ((_distance = vec_ops::distance(_btr->pos(), _com)) > _furthest_distance)
                _furthest_distance = _distance;
            distances.emplace_back(_btr->get_id(), _distance);
        }
        T mid_dist = _furthest_distance/2.0l;
        short b_range1 = _i.b - _h.b;
        short g_range1 = _i.g - _h.g;
        short r_range1 = _i.r - _h.r;
        short b_range2 = _c.b - _i.b;
        short g_range2 = _c.g - _i.g;
        short r_range2 = _c.r - _i.r;
        color _col;
        T _frac;
        std::map<uint64_t, color> _ret_map;
        for (const auto &[_id, _dist] : distances) {
            if (_dist < mid_dist) {
                _frac = _dist/mid_dist;
                _col.b = _h.b + b_range1*_frac;
                _col.g = _h.g + g_range1*_frac;
                _col.r = _h.r + r_range1*_frac;
            } else {
                _frac = _dist/mid_dist - 1;
                _col.b = _i.b + b_range2*_frac;
                _col.g = _i.g + g_range2*_frac;
                _col.r = _i.r + r_range2*_frac;
            }
            _ret_map.emplace(_id, _col);
        }
        return _ret_map;
    }
    template <typename T, typename btrkFuncT, typename funcT, typename ...Args> requires (std::is_fundamental<T>::value)
    std::pair<T, T> run_sim(const char *target_dir,
                            const char *log_path,
                            const std::streamsize &prec,
                            const body<T, T, T, 0>& central_body,
                            long double dt,
                            uint64_t iterations,
                            uint64_t num_reps,
                            const btrkFuncT &btrk_func,
                            const funcT &func,
                            Args ...args) {
        if (!target_dir || !*target_dir)
            throw std::invalid_argument{"Error: target directory cannot be null or empty.\n"};
        // stat directory - delete contents if exists
        if (mkdir(target_dir, S_IRWXU | S_IRWXG | S_IRWXO) == -1)
            throw std::invalid_argument{"Error: directory could not be created.\n"};
        std::ofstream log_file{log_path, std::ios_base::out | std::ios_base::trunc};
        if (!log_file)
            throw std::invalid_argument{"Error: log file could not be opened.\n"};
#ifdef RANDOM_COMET
        auto [_sys, pf, crad] =
#else
        auto [_sys, crad] =
#endif
        func(args...);
        log_file.precision(prec);
        log_file << "Effective radius: " << crad << " m" << "\nBulk density: " <<
        gtd::comet_bd(crad, _sys.num_bodies()*_sys.front().mass()) << " kg/m^3" <<
#ifdef RANDOM_COMET
        "\nPacking fraction: " << pf <<
#endif
        std::endl;
        _sys.add_body(central_body);
        _sys.set_iterations(iterations);
        _sys.set_timestep(dt);
        body_tracker<T, T, T, 0> _btrk{_sys, btrk_func};
        log_file << "Total number of bodies: " << _sys.num_bodies() << "\nNumber of bodies in comet: " <<
        _btrk.num_bods() << std::endl;
        long double base = log10l(num_reps);
        unsigned short num_digits = (unsigned short) ceill(base);
        size_t dir_len = strlen_c(target_dir);
        char *nsys_path;
        char *_seps_path;
        char *_cseps_path;
        if (*(target_dir + dir_len - 1) == '/') {
            nsys_path = new char[14 + num_digits + dir_len]{};
            _seps_path = new char[16 + dir_len]{};
            _cseps_path = new char[17 + dir_len]{};
            strcpy_c(nsys_path, target_dir);
            strcpy_c(_seps_path, target_dir);
            strcpy_c(_cseps_path, target_dir);
        } else {
            std::ptrdiff_t offset = dir_len++;
            nsys_path = new char[14 + num_digits + dir_len]{};
            _seps_path = new char[16 + dir_len]{};
            _cseps_path = new char[17 + dir_len]{};
            strcpy_c(nsys_path, target_dir);
            strcpy_c(_seps_path, target_dir);
            strcpy_c(_cseps_path, target_dir);
            *(nsys_path + offset) = '/';
            *(_seps_path + offset) = '/';
            *(_cseps_path + offset) = '/';
        }
        strcat_c(_seps_path, "separations.sep");
        strcat_c(_cseps_path, "separations.csep");
        strcat_c(nsys_path, "DataFile");
        char *ptr = nsys_path + 8 + dir_len;
        unsigned short _count = 0;
        for (; _count < num_digits; ++_count)
            *ptr++ = '0';
        strcpy_c(ptr, ".nsys");
        unsigned short num_zeros;
        uint64_t to_conv;
        char *last_digit = nsys_path + dir_len + 7 + num_digits;
        uint64_t _i = 0;
        long double time_per_step = iterations*dt;
        vector3D<T> _com;
        vector3D<T> _com_vel;
        std::vector<T> _seps;
        std::vector<T> _cseps;
        goto start_loop;
        for (; _i <= num_reps; ++_i) {
            ptr = last_digit;
            to_conv = _i;
            while (to_conv > 0) {
                *ptr-- = to_conv % 10 + 48;
                to_conv /= 10;
            }
            _sys.evolve(sys::leapfrog_kdk);
            start_loop:
            _sys.to_nsys(nsys_path);
            _com = _btrk.com();
            _com_vel = _btrk.com_vel();
            log_file << "Frame: " << _i << '/' << num_reps << "\nCOM: " << _com << " m" <<
            "\nCOM velocity: " << _com_vel << " m/s\nCOM dist. to c.p.: " <<
            vec_ops::distance(_com, central_body) << " m\n" << "\nCOM speed: " << _com_vel.magnitude() <<
            " m/s\nMean sep.: " << _seps.emplace_back(_btrk.mean_sep()) << " m\nMean dist. to COM: " <<
            _cseps.emplace_back(_btrk.mean_dist_to(_com)) << " m\n"
            "Time: " << time_per_step*_i << " seconds" << std::endl;
        }
        unsigned short T_size = (unsigned short) sizeof(T);
        uint64_t tot_frames = num_reps + 1;
        if (FILE *fp = fopen(_seps_path, "wb"); fp != nullptr) {
            if constexpr (std::floating_point<T>)
                fwrite("FP", sizeof(char), 2, fp);
            else {
                if constexpr (std::signed_integral<T>)
                    fwrite("UI", sizeof(char), 2, fp);
                else
                    fwrite("SI", sizeof(char), 2, fp);
            }
            fwrite(&T_size, sizeof(unsigned short), 1, fp);
            fwrite(&tot_frames, sizeof(uint64_t), 1, fp);
            fwrite(_seps.data(), sizeof(T), tot_frames, fp);
            fclose(fp);
        } else
            log_file << "Error writing mean separation data.\n";
        if (FILE *fp = fopen(_cseps_path, "wb"); fp != nullptr) {
            if constexpr (std::floating_point<T>)
                fwrite("FP", sizeof(char), 2, fp);
            else {
                if constexpr (std::signed_integral<T>)
                    fwrite("UI", sizeof(char), 2, fp);
                else
                    fwrite("SI", sizeof(char), 2, fp);
            }
            fwrite(&T_size, sizeof(unsigned short), 1, fp);
            fwrite(&tot_frames, sizeof(uint64_t), 1, fp);
            fwrite(_cseps.data(), sizeof(T), tot_frames, fp);
            fclose(fp);
        } else
            log_file << "Error writing mean distance to COM data.\n";
        log_file.close();
        return {_seps.back(), _cseps.back()}; // return final mean sep. and final mean dist. to COM of comet
    }
}
#endif
