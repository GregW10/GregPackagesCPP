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
                    const vec3& c_pos, const vec3& c_vel, bool random_comet = false, const char *mass_units = "kg",
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
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, bool prog, bool merge, int coll, uint64_t mF, 
              uint64_t fF, bool binary>
    std::map<uint64_t, color> col_dist(const body_tracker<M, R, T, mF> &btrk, const color &_h, const color &_c) {
        if (!sys.num_bodies())
	    return {};
	vector3D<T> _com = btrk.com();
        T _furthest_distance{};
        T _distance;
        std::vector<std::pair<uint64_t, T>> distances;
	for (const body<M, R, T, mF>* &_btr : btrk) {
            if ((_distance = vec_ops::distance((*_btr)->pos(), _com)) > _furthest_distance)
                _furthest_distance = _distance;
            distances.emplace_back((*_btr)->get_id(), _distance);
        }
        unsigned short b_range = _c.b - _h.b;
        unsigned short g_range = _c.g - _h.g;
        unsigned short r_range = _c.r - _h.r;
        color _col;
        T _frac;
        std::map<uint64_t, color> _ret_map;
        for (const auto &[_id, _dist] : distances) {
            _frac = _dist/_furthest_distance;
            _col.b = _h.b + b_range*(_frac);
            _col.g = _h.g + g_range*(_frac);
            _col.r = _h.r + r_range*(_frac);
            _ret_map.emplace(_id, _col);
        }
        return _ret_map;
    }
}
#endif
