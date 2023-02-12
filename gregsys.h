#ifndef GREGSYS_H
#define GREGSYS_H

#include "gregbod.h"

namespace gtd {
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              bool prog = false, bool mergeOverlappingBodies = false, int collisions = 0,
              ull_t memFreq = 0, ull_t fileFreq = 1>
    class system { // G, below, has not been made static as it varies between objects, depending on units passed
    public:
        static constexpr int two_body = 1; // static constants to determine which integration method to perform
        static constexpr int euler = 2;
        static constexpr int modified_euler = 4;
        static constexpr int midpoint = 8;
        static constexpr int leapfrog_kdk = 16;
        static constexpr int leapfrog_dkd = 32;
        static constexpr int rk4 = 64;
        static constexpr int rk3_8 = 128;
        static constexpr int barnes_hut = 65536; // static constants to determine the force approx. method to use
        static constexpr int fast_multipole = 131072;
        // static inline bool merge_if_overlapped = false;
        static constexpr int no_coll_check = 0;
        static constexpr int overlap_coll_check = 3;
        static constexpr int pred_coll_check = 7;
    private:
        long double G = 66743; // Newtonian constant of Gravitation (* 10^15 m^3 kg^-1 s^-2)
        using bod_t = body<M, R, T, memFreq>;
        using vec_size_t = typename std::vector<bod_t>::size_type;
        using sys_t = system<M, R, T, prog, mergeOverlappingBodies, collisions, memFreq, fileFreq>;
        using mom_t = decltype(M{}*T{});
        std::vector<bod_t> bods; // not using set or map as need fast random access to bodies
        std::set<bod_t, std::less<>> del_bods; // set to store bodies removed from system after collision mergers
        long double dt = 1;
        long double half_dt = dt/2; // I gave it its own variable, since it is a commonly used quantity
        ull_t iterations = 1000;
        long double prev_dt{}; // used in methods called after evolve(), since could be changed by setter
        ull_t prev_iterations{}; // same here
        mom_t min_tot_com_mom{pow(10, 10)}; // minimum sum of magnitudes of COM momenta for two bodies to merge
        // using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        std::map<ull_t, std::vector<T>> pe; // to store potential energies of bodies
        std::map<ull_t, std::vector<T>> energy; // total energy for each body at each iteration (KE + PE)
        //std::map<unsigned long long, std::vector<T>> del_ke; // map to store KE vectors of deleted bodies from mergers
        /* it would have been nice to use std::maps to store potential and total energies for all bodies, as it would
         * have allowed lookup by body_id, but its lookup complexity is O(log(N)), whereas accessing an element in a
         * std::vector by index is O(1) */
        std::vector<T> tot_pe; // total potential energy for the entire system at each iteration
        std::vector<T> tot_ke; // total kinetic energy for the entire system at each iteration
        std::vector<T> tot_e; // total energy for the entire system at each iteration (KE + PE)
        bool evolved = false;
        static inline String def_path = get_home_path<String>() + FILE_SEP + "System_Trajectories_";
        static inline unsigned long long to_ull(String &&str) {
            if (!str.isnumeric())
                return -1;
            // str.strip();
            unsigned long long total = 0;
            for (const char &c : str) {
                total *= 10;
                total += c - 48;
            }
            return total;
        }
        void parse_units_format(String &&str) {
            if (!std::regex_match(str.c_str(), std::regex(R"(^\s*(m|M)\s*\d{1,18}\s*:\s*\d{1,18}\s*,\s*(d|D)\s*\d{1,18}\s*:\s*\d{1,18}\s*,\s*(t|T)\s*\d{1,18}\s*:\s*\d{1,18}\s*$)"))) {
                str.append_front("The units_format string passed, \"").append_back("\", does not match the format "
                                                                                   "required:\n\"Ma:b,Dc:x,Ty:z\", "
                                                                                   "where 'a', 'b', 'c', 'x', 'y', and "
                                                                                   "'z' represent any integer number of"
                                                                                   " 1-18 characters.");
                throw std::invalid_argument(str.c_str());
            }
            str.strip("MmDdTt ");
            size_t colon_index = str.find(':');
            size_t comma_index = str.find(',');
            unsigned long long m_denom = to_ull(str.substr(0, colon_index)); // denominator
            unsigned long long m_num = to_ull(str.substr(colon_index + 1, comma_index)); // numerator
            str.erase_chars(0, comma_index + 1);
            colon_index = str.find(':');
            comma_index = str.find(',');
            unsigned long long d_num = to_ull(str.substr(0, colon_index));
            unsigned long long d_denom = to_ull(str.substr(colon_index + 1, comma_index));
            str.erase_chars(0, comma_index + 1);
            colon_index = str.find(':');
            comma_index = str.find(',');
            unsigned long long t_denom = to_ull(str.substr(0, colon_index));
            unsigned long long t_num = to_ull(str.substr(colon_index + 1, comma_index));
            G *= m_num*d_num*d_num*d_num*t_num*t_num;
            G /= m_denom*d_denom*d_denom*d_denom*t_denom*t_denom;
            G /= 1000000000000000; // correcting for the original G being 10^(15) times larger than it should be
        }
        void clear_bodies() { // makes sure the trajectories of all bodies are deleted
            // for (bod_t &bod : bods)
            //     bod.clear();
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.clear();});
        }
        void clear_bodies(vec_size_t &as_of) {
            // if (as_of >= bods.size())
            //     return;
            // vec_size_t size = bods.size();
            // while (as_of < size)
            //     bods[as_of++].clear();
            std::for_each(bods.begin() + as_of, bods.end(), [this](bod_t &bod){bod.clear();});
        }
        void clear_bodies(vec_size_t &&as_of) {
            clear_bodies(as_of);
        }
        void check_overlap() {
            vec_size_t size;
            bod_t *outer;
            bod_t *inner;
            vec_size_t i = 0;
            vec_size_t j;
            start:
            size = bods.size();
            for (; i < size; ++i) {
                outer = bods.data() + i;
                for (j = i + 1; j < size; ++j) {
                    inner = bods.data() + j;
                    if ((outer->curr_pos - inner->curr_pos).magnitude() < outer->radius + inner->radius) {
                        if constexpr (!mergeOverlappingBodies) {
                            String str = "The bodies with id=";
                            str.append_back(outer->id).append_back(" and id=").append_back(inner->id);
                            str.append_back(" that were added to this system object overlap.\n");
                            throw overlapping_bodies_error(str.c_str());
                        }
                        *outer += *inner; // merges the two overlapping bodies
                        bods.erase(bods.begin() + j); // thus the number of total bodies is reduced by 1
                        goto start;
                    }
                }
            }
        }
        void cumulative_acc_and_pe(bod_t &b1, bod_t &b2) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
            long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
            b1.acc += ((G*b2.mass_)/(r12_cubed_mag))*r12;
            b2.acc -= ((G*b1.mass_)/(r12_cubed_mag))*r12;
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/r12.magnitude();
            b1.pe += pot_energy;
            b2.pe += pot_energy;
        }
        void cumulative_acc(bod_t &b1, bod_t &b2) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
            long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
            b1.acc += ((G*b2.mass_)/(r12_cubed_mag))*r12;
            b2.acc -= ((G*b1.mass_)/(r12_cubed_mag))*r12;
        }
        void cumulative_pe(bod_t &b1, bod_t &b2) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/r12.magnitude();
            b1.pe += pot_energy;
            b2.pe += pot_energy;
        }
        void take_euler_step() {
            for (bod_t &bod : bods) {
                bod.curr_pos += bod.curr_vel*dt;
                bod.curr_vel += bod.acc*dt;
                bod.add_pos_vel_ke();
            }
        }
        void take_modified_euler_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                      &predicted_vals) {
            unsigned long long count = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += half_dt*(bod.curr_vel + std::get<1>(predicted_vals[count]));
                bod.curr_vel += half_dt*(bod.acc + std::get<2>(predicted_vals[count++]));
                bod.add_pos_vel_ke();
            }
        }
        void take_midpoint_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                &predicted_vals) {
            unsigned long long count = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += dt*std::get<1>(predicted_vals[count]);
                bod.curr_vel += dt*std::get<2>(predicted_vals[count++]);
                bod.add_pos_vel_ke();
            }
        }
        void create_energy_vectors() {
            unsigned long long iters_p1 = iterations + 1;
            // pe.resize(bods_size);
            // energy.resize(bods_size);
            pe.clear();
            energy.clear();
            for (bod_t &bod : bods) {
                // pe[i].reserve(iters_p1);
                // energy[i].reserve(iters_p1);
                pe[bod.id].reserve(iters_p1);
                energy[bod.id].reserve(iters_p1);
                bod.positions.reserve(iters_p1);
                bod.velocities.reserve(iters_p1);
                bod.energies.reserve(iters_p1);
            }
            tot_pe.reserve(iters_p1);
            tot_ke.reserve(iters_p1);
            tot_e.reserve(iters_p1);
        }
        void calc_acc_and_e() {
            static T total_pe;
            static T total_ke;
            static unsigned long long outer;
            static unsigned long long inner;
            static vec_size_t num_bods;
            total_pe = T{0};
            total_ke = T{0};
            for (bod_t &bod : bods) {
                bod.acc.make_zero();
                bod.pe = T{0};
            }
            num_bods = bods.size();
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    cumulative_acc_and_pe(ref, bods[inner]);
                }
                ref.pe /= T{2};
                pe[ref.id].push_back(ref.pe);
                energy[ref.id].push_back(ref.curr_ke + ref.pe);
                total_pe += ref.pe;
                total_ke += ref.curr_ke;
            }
            tot_pe.push_back(total_pe);
            tot_ke.push_back(total_ke);
            tot_e.push_back(total_pe + total_ke);
        }
        void leapfrog_kdk_acc_e_and_step() {
            static T total_pe;
            static T total_ke;
            static unsigned long long inner;
            static unsigned long long outer;
            static vec_size_t num_bods;
            total_pe = T{0};
            total_ke = T{0};
            for (bod_t &bod : bods) {
                bod.acc.make_zero();
                bod.pe = T{0};
            }
            num_bods = bods.size();
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    cumulative_acc_and_pe(ref, bods[inner]);
                }
                ref.pe /= T{2};
                /* KICK for half a step */
                ref.curr_vel += half_dt*ref.acc;
                ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
                pe[ref.id].push_back(ref.pe);
                energy[ref.id].push_back(ref.curr_ke + ref.pe);
                total_pe += ref.pe;
                total_ke += ref.curr_ke;
            }
            tot_pe.push_back(total_pe);
            tot_ke.push_back(total_ke);
            tot_e.push_back(total_pe + total_ke);
        }
        void leapfrog_dkd_acc_e_and_step() {
            static T total_pe;
            static T total_ke;
            static unsigned long long inner;
            static unsigned long long outer;
            static vec_size_t num_bods;
            total_pe = T{0};
            total_ke = T{0};
            for (bod_t &bod : bods) {
                bod.acc.make_zero();
                bod.pe = T{0};
            }
            num_bods = bods.size();
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    cumulative_acc(ref, bods[inner]);
                }
                /* KICK for a full step */
                ref.curr_vel += dt*ref.acc;
                /* DRIFT for half a step */
                ref.curr_pos += half_dt*ref.curr_vel;
                /* the reason it is possible to update the outer body's position within the outer loop (without
                 * affecting the synchronisation of the particles) is because, by here, its effect (at its now-previous
                 * position) on all the other particles in the system has been calculated (all subsequent updates of the
                 * accelerations of the other particles no longer depend on the position of the outer body) */
                ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
                total_ke += ref.curr_ke;
            } /* a second loop is required to compute the potential energies based on the updated positions */
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    cumulative_pe(ref, bods[inner]);
                }
                ref.pe /= T{2};
                pe[ref.id].push_back(ref.pe);
                energy[ref.id].push_back(ref.curr_ke + ref.pe);
                total_pe += ref.pe;
            }
            tot_pe.push_back(total_pe);
            tot_ke.push_back(total_ke);
            tot_e.push_back(total_pe + total_ke);
        }
        void calc_energy() {
            vec_size_t size = bods.size();
            for (bod_t &bod : bods)
                bod.pe = T{0};
            T total_pe{0}; // no static variables as this func. is never called inside a loop
            T total_ke{0};
            unsigned long long inner;
            for (unsigned long long outer = 0; outer < size; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < size; ++inner) {
                    // auto &&pot_energy = (G*ref.mass_*bods[inner].mass_)/(ref.curr_pos-bods[inner].curr_pos).magnitude();
                    // ref.pe -= pot_energy;
                    // bods[inner].pe -= pot_energy;
                    cumulative_pe(ref, bods[inner]);
                }
                ref.pe /= T{2};
                pe[ref.id].push_back(ref.pe);
                energy[ref.id].push_back(ref.curr_ke + ref.pe);
                total_pe += ref.pe;
                total_ke += ref.curr_ke;
            }
            tot_pe.push_back(total_pe);
            tot_ke.push_back(total_ke);
            tot_e.push_back(total_pe + total_ke);
        }
        void s_coll() { // "simple" (ahem, ahem) collision detection and evolution
            static unsigned long long inner;
            static unsigned long long outer;
            static vec_size_t num_bods;
            static std::map<long double, std::tuple<bod_t*, decltype(M{}*T{}), decltype(M{}*T{}), vector3D<T>>>
                    overlapping;
            outer = 0;
            num_bods = bods.size();
            bod_t *bod_i;
            bod_t *merging_bod;
            R rad_dist{};
            mom_t max_mom{min_tot_com_mom};
            mom_t curr_mom{};
            long double min_dist = HUGE_VALL; // usually expands to infinity
            mom_t axis_mom_o;
            mom_t axis_mom_i;
            vec_size_t merging_index;
            for (; outer < num_bods; ++outer) {
                merging_bod = nullptr;
                overlapping.clear();
                bod_t &bod_o = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    bod_i = bods.data() + inner;
                    rad_dist = bod_o.radius + bod_i->radius;
                    vector3D<T> &&r12 = bod_i->curr_pos - bod_o.curr_pos;
                    long double &&dist = r12.magnitude(); // DEAL WITH ZERO DISTANCE CASE
                    r12.x /= dist; r12.y /= dist; r12.z /= dist; // more efficient than calling normalise()
                    if (rad_dist > dist) {
                        vector3D<T> &&com_vel = vel_com(bod_o, *bod_i);
                        axis_mom_o = bod_o.momentum()*r12;
                        axis_mom_i = bod_i->momentum()*r12;
                        if ((curr_mom = axis_mom_o - axis_mom_i - (bod_o.mass_ - bod_i->mass_)*(com_vel*r12))>=max_mom){
                            max_mom = curr_mom;
                            merging_bod = bod_i;
                            merging_index = inner;
                            continue;
                        }
                        if (merging_bod == nullptr) {
                            overlapping.emplace(std::piecewise_construct, std::forward_as_tuple(dist),
                                                std::forward_as_tuple(bod_i, axis_mom_o, axis_mom_i, r12));
                        }
                    }
                }
                if (merging_bod != nullptr) {
                    /* I have opted not to use my += overload for adding a body onto another and updating it (this
                     * would avoid performing 2 deletions from the std::vector) because the position, velocity and KE
                     * data for the body before the merger would be mixed with its new "self" after the merger. In
                     * addition, there would then be inconsistency between which bodies end up in the "deleted" bodies
                     * std::set (retaining all their history) and those that remain in the main std::vector. */
                    bod_t &&merged = bod_o + *merging_bod; // create body that is the result of the merger
                    pe.emplace(merged.id, std::vector<T>{}); // add a std::vector to store its potential energies
                    energy.emplace(merged.id, std::vector<T>{}); // same for its total energy (as it stores its own KE)
                    del_bods.emplace(std::move(bod_o)); // move the merged bodies into the deleted bodies std::set
                    del_bods.emplace(std::move(*bod_i));
                    bods.erase(bods.begin() + outer--); // erase the merged bodies from the std::vector
                    bods.erase(bods.begin() + merging_index - 1); // -1 because outer body was just deleted
                    bods.emplace_back(std::move(merged)); // move the new body into the std::vector storing all bodies
                    --num_bods; // there is a net loss of 1 body
                    continue;
                }
                if (overlapping.size()) {
                    T o_vel;
                    T i_vel;
                    T o_minus_i;
                    T new_o_vel;
                    T new_i_vel;
                    long double avg_rest;
                    for (const auto &[_, tup] : overlapping) {
                        o_vel = std::get<1>(tup)/bod_o.mass_; // recalculating is cheaper than adding them to the tuple
                        i_vel = std::get<2>(tup)/bod_i->mass_; // up above
                        o_minus_i = o_vel - i_vel; // v1 - v2
                        avg_rest = (bod_o.rest_c + std::get<0>(tup)->rest_c)/2;
                        if (o_vel > 0 || i_vel < 0) {
                            if (o_vel <= 0 && i_vel < 0)
                                if (o_vel < i_vel) // bodies are already separating
                                    continue; // case for bodies having passed through each other
                            if (o_vel > 0 && i_vel >= 0)
                                if (o_vel < i_vel)
                                    continue; // case for bodies having passed through each other
                            new_o_vel = (std::get<1>(tup) + std::get<2>(tup) -
                                         std::get<0>(tup)->mass_*avg_rest*o_minus_i)/(bod_o.mass_ +
                                                                                      std::get<0>(tup)->mass_);
                            new_i_vel = (std::get<1>(tup) + std::get<2>(tup) +
                                         bod_o.mass_*avg_rest*o_minus_i)/(bod_o.mass_ + std::get<0>(tup)->mass_);
                            bod_o.curr_vel += (new_o_vel - o_vel)*std::get<3>(tup);
                            std::get<0>(tup)->curr_vel += (new_i_vel - i_vel)*std::get<3>(tup);
                        }
                    }
                }
            }
        }
        static inline bool check_option(int option) noexcept {
            unsigned short loword = option & 0x0000ffff;
            unsigned short hiword = option >> 16;
            return !(loword & (loword - 1)) || !loword || loword > rk4 ||
                   !(hiword & (hiword - 1)) || hiword > 2;
        }
        constexpr void check_coll() {
            static_assert(collisions <= 7 && !(collisions & (collisions + 1)),
                          "Invalid collision-checking option.");
        }
        void print_progress(const unsigned long long &step) const noexcept {
#ifndef _WIN32
            printf(CYAN_TXT_START "Iteration " BLUE_TXT_START "%llu" RED_TXT_START "/" MAGENTA_TXT_START
                   "%llu\r", step, iterations);
#else
            printf("Iteration %llu/%llu\r", step, iterations);
#endif
        }
        static inline void print_conclusion(const char *method, const time_t &start_time) {
            time_t total = time(nullptr) - start_time;
#ifndef _WIN32
            std::cout << RESET_TXT_FLAGS;
            printf(BLACK_TXT("\n--------------------Done--------------------\n")
                   UNDERLINED_TXT(GREEN_TXT("%s"))WHITE_TXT(" - time elapsed: ")
                   BOLD_TXT(YELLOW_TXT("%zu"))WHITE_TXT(" second%c\n"), method, total, "s"[total == 1]);
#else
            printf("\n--------------------Done--------------------\n"
                               "Midpoint method - time elapsed: %zu second%c\n", total, "s"[total == 1]);
#endif
        }
    public:
        system(const std::initializer_list<bod_t> &list) : bods{list} {
            check_overlap();
            check_coll();
            parse_units_format("M1:1,D1:1,T1:1");
            clear_bodies();
        }
        system(std::initializer_list<bod_t> &&list) : bods{std::move(list)} {
            check_overlap();
            check_coll();
            parse_units_format("M1:1,D1:1,T1:1");
            clear_bodies();
        }
        explicit system(long double timestep = 1, unsigned long long num_iterations = 1000,
                        const char *units_format = "M1:1,D1:1,T1:1") :
                dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            check_coll();
            parse_units_format(units_format);
        }
        system(const std::vector<body<M, R, T>> &bodies, long double timestep = 1,
               unsigned long long num_iterations = 1000, const char *units_format = "M1:1,D1:1,T1:1") :
                bods{bodies}, dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            check_overlap();
            check_coll();
            parse_units_format(units_format);
            clear_bodies();
        }
        system(std::vector<body<M, R, T>> &&bodies, long double timestep = 1,
               unsigned long long num_iterations = 1000, const char *units_format = "M1:1,D1:1,T1:1") :
                bods{std::move(bodies)}, dt{timestep}, iterations{num_iterations} {
            check_overlap();
            check_coll();
            parse_units_format(units_format);
            clear_bodies();
        }
        /* copy constructors are made to only copy the variables seen below: it does not make sense for a new system
         * object (even when copy-constructed) to be "evolved" if it has not gone through the evolution itself */
        template <bool prg, bool mrg, int coll> // so a system can be constructed from another without same checks
        system(const system<M, R, T, prg, mrg, coll> &other) :
                bods{other.bods}, dt{other.dt}, iterations{other.iterations}, G{other.G} {
            check_overlap();
            check_coll();
            clear_bodies();
        }
        template <bool prg, bool mrg, int coll>
        system(system<M, R, T, prg, mrg, coll> &&other) :
                bods{std::move(other.bods)}, dt{other.dt}, iterations{other.iterations}, G{other.G} {
            check_overlap();
            check_coll();
            clear_bodies();
        }
        vec_size_t num_bodies() {
            return bods.size();
        }
        bool set_iterations(const unsigned long long &number) noexcept {
            if (!number) // cannot have zero iterations
                return false;
            iterations = number;
            return true;
        }
        bool set_iterations(const unsigned long long &&number) noexcept {
            return set_iterations(number);
        }
        bool set_timestep(const long double &delta_t) noexcept {
            if (delta_t <= 0)
                return false;
            dt = delta_t;
            half_dt = dt/2;
            return true;
        }
        bool set_timestep(const long double &&delta_t) noexcept {
            return set_timestep(delta_t);
        }
        sys_t &add_body(const bod_t &bod) {
            bods.push_back(bod);
            bods.back().clear(); // all bodies within a system object must start out without an evolution
            check_overlap();
            return *this;
        }
        sys_t &add_body(bod_t &&bod) {
            bod.clear();
            bods.push_back(std::move(bod));
            check_overlap();
            return *this;
        }
        template <typename ...Args>
        sys_t &emplace_body(Args&& ...args) { // to allow a body to be constructed in-place, within the system object
            bods.emplace_back(std::forward<Args>(args)...);
            return *this;
        }
        sys_t &add_bodies(const std::vector<bod_t> &bodies) {
            if (!bodies.size())
                return *this;
            vec_size_t index = bods.size();
            bods.insert(bods.end(), bodies.begin(), bodies.end());
            clear_bodies(index);
            check_overlap();
            return *this;
        }
        sys_t &add_bodies(std::vector<bod_t> &&bodies) {
            if (!bodies.size())
                return *this;
            for (auto &b : bodies) {
                b.clear();
                bods.emplace_back(std::move(b));
            }
            check_overlap();
            return *this;
        }
        const bod_t &get_body(unsigned long long id) {
            for (const auto &b : *this)
                if (b.id == id)
                    return b;
            throw std::invalid_argument("The id passed does not correspond to any body present in this system "
                                        "object.\n");
        }
        bool remove_body(unsigned long long id) {
            auto end_it = bods.cend();
            for (typename std::vector<bod_t>::const_iterator it = bods.cbegin(); it < end_it; ++it) {
                if ((*it).id == id) {
                    bods.erase(it);
                    return true;
                }
            }
            return false;
        }
        void reset() {
            // for (bod_t &bod : bods)
            //     bod.reset(true);
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.reset(true);});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
        }
        /* reset() deletes the evolution of the bodies AND returns them to their original positions, velocities, and
         * kinetic energies, whilst clear_evolution() deletes the evolution BUT retains the current positions,
         * velocities and kinetic energies of the bodies. */
        void clear_evolution() {
            // for (bod_t &bod : bods)
            //     bod.clear();
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.clear();});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
        }
        bool evolve(int integration_method = leapfrog_kdk) {
            if (!bods.size() || !check_option(integration_method))
                return false;
            time_t start = time(nullptr);
            unsigned long long steps = 0;
            if (evolved)
                this->clear_evolution();
            vec_size_t num_bods = bods.size();
            this->create_energy_vectors();
#ifndef _WIN32
            if constexpr (prog)
                std::cout << BOLD_TXT_START;
#endif
            if ((integration_method & two_body) == two_body) {
                if (bods.size() != 2)
                    throw two_body_error();
                /* no force approximation techniques performed here, as the integration is just for two bodies */
            }
            else if ((integration_method & euler) == euler) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    while (steps++ < iterations) {
                        if constexpr (collisions == 7) {

                        }
                        this->calc_acc_and_e();
                        this->take_euler_step();
                        /* by having "prog" as a template parameter, it avoids having to re-evaluate whether it is true
                         * or not within the loop, since the "if constexpr ()" branches get rejected at compile-time */
                        if constexpr (collisions == 3)
                            this->s_coll();
                        if constexpr (prog)
                            this->print_progress(steps);
                    }
                    if constexpr (prog)
                        this->print_conclusion("Euler", start);
                }
            }
            else if ((integration_method & modified_euler) == modified_euler) {
                /* The below std::vector will hold the euler-predicted values for position, velocity and acceleration,
                 * respectively. These are subsequently used to determine the corrected values. */
                std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>> predicted{num_bods};
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    unsigned long long outer;
                    unsigned long long inner;
                    while (steps++ < iterations) {
                        this->calc_acc_and_e();
                        for (outer = 0; outer < num_bods; ++outer) {
                            bod_t &ref = bods[outer];
                            std::get<0>(predicted[outer]) = ref.curr_pos + dt*ref.curr_vel;
                            std::get<1>(predicted[outer]) = ref.curr_vel + dt*ref.acc;
                        }
                        for (outer = 0; outer < num_bods; ++outer) {
                            bod_t &ref = bods[outer];
                            std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &pred_outer = predicted[outer];
                            for (inner = outer + 1; inner < num_bods; ++inner) {
                                bod_t &ref_i = bods[inner];
                                vector3D<T> &&r12 = std::get<0>(predicted[inner]) - std::get<0>(pred_outer);
                                long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
                                std::get<2>(pred_outer) += ((G*ref_i.mass_)/(r12_cubed_mag))*r12;
                                std::get<2>(predicted[inner]) -= ((G*ref.mass_)/(r12_cubed_mag))*r12;
                            }
                        }
                        this->take_modified_euler_step(predicted);
                        for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup : predicted)
                            std::get<2>(tup).make_zero();
                        if constexpr (prog)
                            this->print_progress(steps);
                    }
                    if constexpr (prog)
                        this->print_conclusion("Modified Euler", start);
                }
            }
            else if ((integration_method & midpoint) == midpoint) {
                std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>> predicted{num_bods};
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    unsigned long long outer;
                    unsigned long long inner;
                    while (steps++ < iterations) {
                        this->calc_acc_and_e();
                        for (outer = 0; outer < num_bods; ++outer) {
                            bod_t &ref = bods[outer];
                            std::get<0>(predicted[outer]) = ref.curr_pos + half_dt*ref.curr_vel;
                            std::get<1>(predicted[outer]) = ref.curr_vel + half_dt*ref.acc;
                        }
                        for (outer = 0; outer < num_bods; ++outer) {
                            bod_t &ref = bods[outer];
                            std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &pred_outer = predicted[outer];
                            for (inner = outer + 1; inner < num_bods; ++inner) {
                                bod_t &ref_i = bods[inner];
                                vector3D<T> &&r12 = std::get<0>(predicted[inner]) - std::get<0>(pred_outer);
                                long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
                                std::get<2>(pred_outer) += ((G*ref_i.mass_)/(r12_cubed_mag))*r12;
                                std::get<2>(predicted[inner]) -= ((G*ref.mass_)/(r12_cubed_mag))*r12;
                            }
                        }
                        this->take_midpoint_step(predicted);
                        for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup : predicted)
                            std::get<2>(tup).make_zero();
                        if constexpr (prog)
                            this->print_progress(steps);
                    }
                    if constexpr (prog)
                        this->print_conclusion("Midpoint method", start);
                }
            }
            else if ((integration_method & leapfrog_kdk) == leapfrog_kdk) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    this->calc_acc_and_e(); // need an initial acceleration
                    while (steps++ < iterations) {
                        for (bod_t &bod : bods) {
                            /* KICK for half a step */
                            bod.curr_vel += half_dt*bod.acc;
                            /* DRIFT for a full step */
                            bod.curr_pos += dt*bod.curr_vel;
                        }
                        /* update accelerations and energies based on new positions and KICK for half a step */
                        this->leapfrog_kdk_acc_e_and_step();
                        if constexpr (collisions == 3)
                            this->s_coll();
                        if constexpr (prog)
                            this->print_progress(steps);
                    }
                    if constexpr (prog)
                        this->print_conclusion("Leapfrog KDK", start);
                }
                goto bye;
            }
            else if ((integration_method & leapfrog_dkd) == leapfrog_dkd) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    this->calc_energy();
                    while (steps++ < iterations) {
                        for (bod_t &bod : bods)
                            /* DRIFT for half a step */
                            bod.curr_pos += half_dt*bod.curr_vel;
                        /* update accelerations, then KICK for a full step and DRIFT for half a step, then update E */
                        this->leapfrog_dkd_acc_e_and_step();
                        if constexpr (prog)
                            this->print_progress(steps);
                    }
                    if constexpr (prog)
                        this->print_conclusion("Leapfrog DKD", start);
                }
                goto bye;
            }
            else if ((integration_method & rk4) == rk4) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {

                }
            }
            else if ((integration_method & rk3_8) == rk3_8) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {

                }
            }
            this->calc_energy();
            bye:
            evolved = true;
            prev_dt = dt;
            prev_iterations = iterations;
            return true;
        }
        std::pair<T*, long double> *pe_extrema() const {
            /* This method returns a pointer to a 2-element array of std::pairs. The first pair concerns the minima and
             * the second pair concerns the maxima. The first element of the pairs is a pointer to a 2-element array of
             * type T. The first element of the array is the actual extreme value, and the second is the difference
             * between the extreme value and starting PE. The second element of the pairs is the time at which this
             * extreme occurs. */
            if (bods.empty())
                throw empty_system_error("pe_diff() cannot be called on empty system objects (no bodies present).\n");
            if (!evolved)
                throw no_evolution_error("pe_diff() can only be called after a system object has been evolved.\n");
            static T minima[2]; // first the min. PE, then the difference between it and the starting PE
            static T maxima[2]; // same but for max.
            static std::pair<T*, long double> min_max[2] = {{minima, 0}, {maxima, 0}};
            unsigned long long min_index{};
            unsigned long long max_index{};
            unsigned long long count{};
            const T &first_element = tot_pe[0]; // avoids 3 function calls
            T min = first_element;
            T max = first_element;
            for (const T &pot_e : tot_pe) {
                if (pot_e < min) {
                    min = pot_e;
                    min_index = count;
                }
                if (pot_e > max) {
                    max = pot_e;
                    max_index = count;
                }
                ++count;
            }
            *minima = min; // assign extrema
            *maxima = max;
            *(minima + 1) = first_element - min; // assign differences from initial PE
            *(maxima + 1) = max - first_element;
            min_max->second = min_index*prev_dt; // assign times at which extrema occurred
            (min_max + 1)->second = max_index*prev_dt;
            return min_max;
        }
        std::pair<T*, long double> *ke_extrema() const {
            /* Same as the above method, but for KE. */
            if (bods.empty())
                throw empty_system_error("pe_diff() cannot be called on empty system objects (no bodies present).\n");
            if (!evolved)
                throw no_evolution_error("pe_diff() can only be called after a system object has been evolved.\n");
            static T minima[2]; // first the min. PE, then the difference between it and the starting PE
            static T maxima[2]; // same but for max.
            static std::pair<T*, long double> min_max[2] = {{minima, 0}, {maxima, 0}};
            unsigned long long min_index{};
            unsigned long long max_index{};
            unsigned long long count{};
            const T &first_element = tot_ke[0];
            T min = first_element;
            T max = first_element;
            for (const T &kin_e : tot_ke) {
                if (kin_e < min) {
                    min = kin_e;
                    min_index = count;
                }
                if (kin_e > max) {
                    max = kin_e;
                    max_index = count;
                }
                ++count;
            }
            *minima = min;
            *maxima = max;
            *(minima + 1) = first_element - min;
            *(maxima + 1) = max - first_element;
            min_max->second = min_index*prev_dt;
            (min_max + 1)->second = max_index*prev_dt;
            return min_max;
        }
        std::pair<T*, long double> *tot_e_extrema() const {
            /* Same as the above method, but for KE. */
            if (bods.empty())
                throw empty_system_error("pe_diff() cannot be called on empty system objects (no bodies present).\n");
            if (!evolved)
                throw no_evolution_error("pe_diff() can only be called after a system object has been evolved.\n");
            static T minima[2]; // first the min. PE, then the difference between it and the starting PE
            static T maxima[2]; // same but for max.
            static std::pair<T*, long double> min_max[2] = {{minima, 0}, {maxima, 0}};
            unsigned long long min_index{};
            unsigned long long max_index{};
            unsigned long long count{};
            const T &first_element = tot_e[0];
            T min = first_element;
            T max = first_element;
            for (const T &_e : tot_e) {
                if (_e < min) {
                    min = _e;
                    min_index = count;
                }
                if (_e > max) {
                    max = _e;
                    max_index = count;
                }
                ++count;
            }
            *minima = min;
            *maxima = max;
            *(minima + 1) = first_element - min;
            *(maxima + 1) = max - first_element;
            min_max->second = min_index*prev_dt;
            (min_max + 1)->second = max_index*prev_dt;
            return min_max;
        }
        auto begin() {
            return bods.begin();
        }
        auto end() {
            return bods.end();
        }
        auto begin() const {
            return bods.cbegin();
        }
        auto end() const {
            return bods.cend();
        }
        auto cbegin() const {
            return bods.cbegin();
        }
        auto cend() const {
            return bods.cend();
        }
        auto rbegin() {
            return bods.rbegin();
        }
        auto rend() {
            return bods.rend();
        }
        auto rbegin() const {
            return bods.crbegin();
        }
        auto rend() const {
            return bods.crend();
        }
        auto crbegin() const {
            return bods.crbegin();
        }
        auto crend() const {
            return bods.crend();
        }
        bool write_trajectories(const String &path = def_path, bool binary = false) {
            /* Writes the trajectories (historically) of all the bodies to a file, including all the energy values. For
             * text files, units are not written as these will depend entirely on the units_format string passed to the
             * system<M, R, T> object in its constructor. */
            if (&path == &def_path)
                def_path.append_back(get_date_and_time()).append_back(".csv");
            if (!evolved || bods.size() == 0)
                return false;
            if (binary) {} // needs work
            std::ofstream out(path.c_str(), std::ios_base::trunc);
            if (!out.good())
                return false;
            unsigned long long count = 0;
            unsigned long long i;
            for (const bod_t &bod : bods) {
                out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                    << "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy,"
                       "potential_energy,total_energy\r\n";
                for (i = 0; i <= prev_iterations; ++i)
                    out << i*prev_dt << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                        << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                        << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                        << energy[bod.id][i] << "\r\n";
                out << ",,,,,,,,,\r\n";
                ++count;
            }
            out << "time_elapsed,system_PE,system_KE,system_E\r\n";
            for (i = 0; i <= prev_iterations; ++i)
                out << i*prev_dt << ',' << tot_pe[i] << ',' << tot_ke[i] << ',' << tot_e[i] << "\r\n";
            out << ",,,,,,,,,\r\n";
            out << "energy_type,min,max,min_abs_diff_from_start_val,"
                   "max_abs_diff_from_start_val,time_of_min,time_of_max\r\n";
            std::pair<T*, long double> *ptr = this->pe_extrema();
            out << "system_PE," << *ptr->first << ',' << *((ptr + 1)->first) << ',' << *(ptr->first + 1) << ','
                << *((ptr + 1)->first + 1) << ',' << ptr->second << ',' << (ptr + 1)->second << "\r\n";
            ptr = this->ke_extrema();
            out << "system_KE," << *ptr->first << ',' << *((ptr + 1)->first) << ',' << *(ptr->first + 1) << ','
                << *((ptr + 1)->first + 1) << ',' << ptr->second << ',' << (ptr + 1)->second << "\r\n";
            ptr = this->tot_e_extrema();
            out << "system_E," << *ptr->first << ',' << *((ptr + 1)->first) << ',' << *(ptr->first + 1) << ','
                << *((ptr + 1)->first + 1) << ',' << ptr->second << ',' << (ptr + 1)->second << "\r\n";
            out.close();
            if (&path == &def_path)
                def_path.erase_chars(def_path.get_length() - 29);
            // std::cout << "potential energy size: " << pe.size() << std::endl;
            // std::cout << "energy size: " << energy.size() << std::endl;
            // std::cout << "potential energy v size: " << pe[0].size() << std::endl;
            // std::cout << "energy v size: " << energy[0].size() << std::endl;
            // std::cout << "Total potential energy size: " << tot_pe.size() << std::endl;
            // std::cout << "Total kinetic energy size: " << tot_ke.size() << std::endl;
            // std::cout << "Total energy size: " << tot_e.size() << std::endl;
            return true;
        }
        const bod_t &operator[](vec_size_t index) {
            if (index >= bods.size())
                throw std::out_of_range("The requested body does not exist (index out of range).\n");
            return bods[index];
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg, bool mrg, int coll, ull_t mF, ull_t fF>
        sys_t &operator=(const system<m, r, t, prg, mrg, coll, mF, fF> &other) {
            if (this == &other)
                return *this;
            this->dt = other.dt;
            this->iterations = other.iterations;
            this->G = other.G;
            this->bods = other.bods;
            this->clear_evolution();
            return *this;
        }
        sys_t &operator=(sys_t &&other) noexcept {
            if (this == &other)
                return *this;
            this->dt = std::move(other.dt);
            this->iterations = std::move(other.iterations);
            this->G = std::move(other.G);
            this->bods = std::move(other.bods);
            this->clear_evolution();
            return *this;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg, bool mrg, int coll, ull_t mF, ull_t fF>
        friend std::ostream &operator<<(std::ostream&, const system<m, r, t, prg, mrg, coll, mF, fF>&);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg1, bool prg2, bool mrg1, bool mrg2,
                  int c1, int c2, ull_t mF1, ull_t mF2, ull_t fF1, ull_t fF2>
        friend system<m, r, t, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2)>
        operator+(const system<m, r, t, prg1, mrg1, c1, mF1, fF1>&,
                  const system<m, r, t, prg2, mrg2, c2, mF2, fF2>&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, ull_t, ull_t>
        friend class system;
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                isNumWrapper, bool>
        friend class astro_scene;
    }; // all these template parameters are driving me coocoo
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, bool prg, bool mrg, int coll, ull_t mF, ull_t fF>
    std::ostream &operator<<(std::ostream &os, const system<M, R, T, prg, mrg, coll, mF, fF> &sys) {
        typename std::vector<body<M, R, T, mF>>::size_type count = sys.bods.size();
        os << "[gtd::system@" << &sys << ",num_bodies=" << count;
        if (!count)
            return os << ']';
        os << ",bodies:";
        count = 0;
        for (const body<M, R, T, mF> &b : sys)
            os << "\n body_" << count++ << '=' << b;
        return os << ']';
    }
    template<isNumWrapper M, isNumWrapper R, isNumWrapper T, bool prg1, bool prg2, bool mrg1, bool mrg2, int c1, int c2,
            ull_t mF1, ull_t mF2, ull_t fF1, ull_t fF2>
    system<M, R, T, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2)>
    operator+(const system<M, R, T, prg1, mrg1, c1, mF1, fF1> &sys1,
              const system<M, R, T, prg2, mrg2, c2, mF2, fF2> &sys2) {
        typename std::vector<body<M, R, T, mF1>>::size_type size1 = sys1.bods.size();
        typename std::vector<body<M, R, T, mF2>>::size_type size2 = sys2.bods.size();
        if ((!size1 && !size2) || sys1.G != sys2.G)
            return {{}};
        if (!size1)
            return {sys2};
        if (!size2)
            return {sys1};
        system<M, R, T, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2)>
                ret_sys(sys1.bods, (sys1.dt + sys2.dt)/2.0l, (sys1.iterations + sys2.iterations)/2);
        // for (typename std::vector<body<M, R, T>>::size_type i = 0; i < size1; ++i) {
        //     for (typename std::vector<body<M, R, T>>::size_type j = 0; j < size2; ++j) {
        //         if ((sys1.bods[i].curr_pos-sys2.bods[j].curr_pos).magnitude()<sys1.bods[i].radius+sys2.bods[j].radius) {
        //             if (!system<M, R, T>::merge_if_overlapped)
        //                 throw overlapping_bodies_error();
        //             ret_sys.bods[i];
        //         }
        //     }
        // }
        ret_sys.G = sys1.G;
        ret_sys.add_bodies(sys2.bods);
        return ret_sys;
    }
    // std::set<unsigned long long> body_counter::ids;
    // unsigned long long body_counter::count = 0;
    typedef system<long double, long double, long double> sys;
    typedef system<long double, long double, long double, true> sys_p;
    typedef system<long double, long double, long double, false, true> sys_m;
    typedef system<long double, long double, long double, false, false, 3> sys_c;
    typedef system<long double, long double, long double, false, false, 7> sys_C;
    typedef system<long double, long double, long double, true, true> sys_pm;
    typedef system<long double, long double, long double, true, false, 3> sys_pc;
    typedef system<long double, long double, long double, true, false, 7> sys_pC;
    typedef system<long double, long double, long double, false, true, 3> sys_mc;
    typedef system<long double, long double, long double, false, true, 7> sys_mC;
    typedef system<long double, long double, long double, true, true, 3> sys_pmc;
    typedef system<long double, long double, long double, true, true, 7> sys_pmC;
    typedef system<long double, long double, long double, false, false, 0, 1, 0> sys_M;
    typedef system<long double, long double, long double, true, false, 0, 1, 0> sys_p_M;
    typedef system<long double, long double, long double, false, true, 0, 1, 0> sys_m_M;
    typedef system<long double, long double, long double, false, false, 3, 1, 0> sys_c_M;
    typedef system<long double, long double, long double, false, false, 7, 1, 0> sys_C_M;
    typedef system<long double, long double, long double, true, true, 0, 1, 0> sys_pm_M;
    typedef system<long double, long double, long double, true, false, 3, 1, 0> sys_pc_M;
    typedef system<long double, long double, long double, true, false, 7, 1, 0> sys_pC_M;
    typedef system<long double, long double, long double, false, true, 3, 1, 0> sys_mc_M;
    typedef system<long double, long double, long double, false, true, 7, 1, 0> sys_mC_M;
    typedef system<long double, long double, long double, true, true, 3, 1, 0> sys_pmc_M;
    typedef system<long double, long double, long double, true, true, 7, 1, 0> sys_pmC_M;
    typedef system<long double, long double, long double, false, false, 0, 1, 1> sys_MF;
    typedef system<long double, long double, long double, true, false, 0, 1, 1> sys_p_MF;
    typedef system<long double, long double, long double, false, true, 0, 1, 1> sys_m_MF;
    typedef system<long double, long double, long double, false, false, 3, 1, 1> sys_c_MF;
    typedef system<long double, long double, long double, false, false, 7, 1, 1> sys_C_MF;
    typedef system<long double, long double, long double, true, true, 0, 1, 1> sys_pm_MF;
    typedef system<long double, long double, long double, true, false, 3, 1, 1> sys_pc_MF;
    typedef system<long double, long double, long double, true, false, 7, 1, 1> sys_pC_MF;
    typedef system<long double, long double, long double, false, true, 3, 1, 1> sys_mc_MF;
    typedef system<long double, long double, long double, false, true, 7, 1, 1> sys_mC_MF;
    typedef system<long double, long double, long double, true, true, 3, 1, 1> sys_pmc_MF;
    typedef system<long double, long double, long double, true, true, 7, 1, 1> sys_pmC_MF;
}
#endif
