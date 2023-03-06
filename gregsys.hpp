#ifndef GREGSYS_H
#define GREGSYS_H

#ifndef __cplusplus
#error "The gregsys.hpp header file is a C++ header file only."
#endif

#include "gregbod.hpp"

#define EMPTY

#define FUNC_TEMPL_SELECT(func, cf, ...) \
if constexpr (collisions == overlap_coll_check || !(memFreq && fileFreq)) { \
    func<false, false>(__VA_ARGS__);            \
    if constexpr (collisions == overlap_coll_check) {                       \
        if constexpr (memFreq && fileFreq) {                           \
            if (!(steps % memFreq) && !(steps % fileFreq)) { \
                this->s_coll<true, true>(); \
            } \
            else if (!(steps % memFreq)) { \
                this->s_coll<true, false>(); \
            } \
            else if (!(steps % fileFreq)) { \
                this->s_coll<false, true>(); \
                cf \
            } \
            else { \
                this->s_coll<false, false>(); \
                cf \
            } \
        }                                \
        else if constexpr (memFreq) { \
            if (steps % memFreq) { \
                this->s_coll<false, false>(); \
                cf \
            } \
            else { \
                this->s_coll<true, false>(); \
            } \
        } \
        else if constexpr (fileFreq) { \
            if (steps % fileFreq) { \
                this->s_coll<false, false>(); \
                cf \
            } \
            this->s_coll<false, true>(); \
            cf \
        } \
    }                                    \
    else { \
        cf                               \
    }\
} \
else { \
    if constexpr (memFreq && fileFreq) { \
        if (!(steps % memFreq) && !(steps % fileFreq)) { \
            func<true, true>(__VA_ARGS__); \
        } \
        else if (!(steps % memFreq)) { \
            func<true, false>(__VA_ARGS__); \
        } \
        else if (!(steps % fileFreq)) { \
            func<false, true>(__VA_ARGS__); \
            cf \
        } \
        else { \
            func<false, false>(__VA_ARGS__); \
            cf \
        } \
    } \
    else if constexpr (memFreq) { \
        if (steps % memFreq) { \
            func<false, false>(__VA_ARGS__); \
            cf \
        } \
        else { \
            func<true, false>(__VA_ARGS__); \
        } \
    } \
    else if constexpr (fileFreq) { \
        if (steps % fileFreq) { \
            func<false, false>(__VA_ARGS__); \
            cf \
        } \
        func<false, true>(__VA_ARGS__); \
        cf \
    } \
}

#define MEM_LOOP \
if constexpr (memFreq) { \
    if (steps % memFreq) \
        for (bod_t &bod : bods) \
            bod.acc.make_zero(); \
    else {       \
        total_pe = total_ke =  T{0}; \
        for (bod_t &bod : bods) { \
            bod.acc.make_zero(); \
            bod.pe = T{0}; \
        } \
    }\
} \
else \
    for (bod_t &bod : bods) \
        bod.acc.make_zero();

namespace gtd {
    class nsys_load_error : public nbody_error {
    public:
        nsys_load_error() : nbody_error{"Error loading N-body data from .nsys file.\n"} {}
        explicit nsys_load_error(const char *msg) : nbody_error{msg} {}
    };
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              bool prog = false, bool mergeOverlappingBodies = false, int collisions = 0,
              uint64_t memFreq = 0, uint64_t fileFreq = 1, bool binaryFile = true>
    class system {
        /* The system class represents a system in 3D space that is composed of bodies which interact with each other
         * gravitationally (and through contact force in case of collisions - if this option is specified using the
         * collisions template parameter (3 or 7)). The bodies in the system are represented by gtd::body objects. */
        /* A system object is capable of evolving, in time, the positions, velocities & energies of the bodies via
         * various different integration methods, each specified by an integral constant (just below). */
        /* The 'M', 'R' & 'T' template parameters are the data types used to represent the bodies' masses, radii, &
         * positions and velocities, respectively. The other (non-type) template parameters are present in order to
         * avoid a small runtime overhead in checking which have been set and/or to avoid duplication of code. The
         * 'prog' parameter determines whether progress is printed to standard output during the main evolution of the
         * bodies. 'mergeOverlappingBodies' determines whether bodies that are added to a system object are merged
         * together if they overlap, or whether an exception is thrown instead. 'collisions' determines how collisions
         * should be predicted and how they should be dealt with (0 = not at all, 3 = simple - based on overlap, 7 =
         * advanced (predictive) collision checking). 'memFreq' represents the frequency with which position, velocity
         * and energy data values should be stored in memory (on the heap) and, similarly, fileFreq represents the
         * frequency with which they should be written to a file. Finally, 'binaryFile' represents whether the data
         * being written to a file is in binary format or text (.csv) format (if fileFreq is non-zero). */
    public:
        static constexpr int two_body = 1; // static constants to determine which integration method to perform
        static constexpr int euler = 2;
        static constexpr int modified_euler = 4;
        static constexpr int midpoint = 8;
        static constexpr int leapfrog_kdk = 16;
        static constexpr int leapfrog_dkd = 32;
        static constexpr int rk4 = 64;
        static constexpr int rk3_8 = 128;
        static constexpr int barnes_hut = 65'536; // static constants to determine the force approx. method to use
        static constexpr int fast_multipole = 131'072;
        static constexpr int no_coll_check = 0;
        static constexpr int overlap_coll_check = 3;
        static constexpr int pred_coll_check = 7;
        static constexpr long double G_SI = 0.000'000'000'066'743; // G in SI units (m^3 kg^-1 s^-2)
    private: // G, below, has not been made a static constant as it varies between objects, depending on units passed
        long double G = 66'743; // Newtonian constant of Gravitation (10^(-15) m^3 kg^-1 s^-2)
        using bod_t = body<M, R, T, memFreq>;
        using vec_size_t = typename std::vector<bod_t>::size_type;
        using sys_t = system<M, R, T, prog, mergeOverlappingBodies, collisions, memFreq, fileFreq, binaryFile>;
        using mom_t = decltype(M{}*T{});
        std::vector<bod_t> bods; // not using set, map or list as need fast random access to bodies
        // std::set<bod_t, std::less<>> del_bods; // set to store bodies removed from system after collision mergers
        /* Here I declare a std::set to store bodies that are removed from the evolution in the case of a merger. I have
         * declared the type of the comparator used by the std::set object to be that of an anonymous lambda type, which
         * will cause the std::pair objects to be sorted first by the iteration at which they were "destroyed", and then
         * by their IDs (for when times of destruction are equal). */
        // std::set<std::pair<ull_t, bod_t>,
        //         decltype([](const std::pair<ull_t, bod_t> &p1, const std::pair<ull_t, bod_t> &p2){
        //             return p1.first == p2.first ? p1.second < p2.second : p1.first < p2.first;
        //         })> del_bods;
        static inline auto pair_bods_func = [](const std::pair<uint64_t, bod_t> &p1,
                                               const std::pair<uint64_t, bod_t> &p2) {
            return p1.first == p2.first ? p1.second < p2.second : p1.first < p2.first;
        };
        union {
            std::set<bod_t, std::less<>> *del_bods_n;
            std::set<std::pair<uint64_t, bod_t>, decltype(pair_bods_func)> *del_bods_m;
        };
        long double dt = 1;
        long double half_dt = dt/2; // I gave it its own variable, since it is a commonly used quantity
        uint64_t iterations = 1'000;
        long double prev_dt{}; // used in methods called after evolve(), since could be changed by setter
        long double time_elapsed{};
        uint64_t prev_iterations{}; // same here
        mom_t min_tot_com_mom{BILLION*10}; // minimum sum of magnitudes of COM momenta for two bodies to merge
        // using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        uint64_t steps{}; // defined as an instance variable since it's required in numerous functions
        T total_pe{}; // these 5 variables are defined as instance variables to avoid their redefinition in many funcs
        T total_ke{};
        uint64_t inner{};
        uint64_t outer{};
        mutable vec_size_t num_bods{};
        std::map<uint64_t, std::vector<T>> pe; // to store potential energies of bodies
        std::map<uint64_t, std::vector<T>> energy; // total energy for each body at each iteration (KE + PE)
        //std::map<unsigned long long, std::vector<T>> del_ke; // map to store KE vectors of deleted bodies from mergers
        /* it would have been nice to use std::maps to store potential and total energies for all bodies, as it would
         * have allowed lookup by body_id, but its lookup complexity is O(log(N)), whereas accessing an element in a
         * std::vector by index is O(1) */
        std::vector<T> tot_pe; // total potential energy for the entire system at each iteration
        std::vector<T> tot_ke; // total kinetic energy for the entire system at each iteration
        std::vector<T> tot_e; // total energy for the entire system at each iteration (KE + PE)
        bool evolved = false; // indicates whether a gtd::system object is in an evolved state
        uint64_t *G_vals = new uint64_t[6];
        union {
            mutable std::ofstream *ostream = nullptr; // to write the system's data to a file
            mutable std::ifstream *istream; // to read data from .nsys files
        };
        static inline String def_path = get_home_path<String>() + FILE_SEP + "System_Trajectories_";
        static inline uint64_t to_ull(String &&str) {
            if (!str.isnumeric())
                return -1;
            uint64_t total = 0;
            for (const char &c : str) {
                total *= 10;
                total += c - 48;
            }
            return total;
        }
        void parse_units_format(String &&str) {
            if (!std::regex_match(str.c_str(), std::regex(R"(^\s*(m|M)\s*(?=\d{1,18}\s*:)\d*[1-9]+\d*\s*:\s*(?=\d{1,18}\s*,)\d*[1-9]+\d*\s*,\s*(d|D)\s*(?=\d{1,18}\s*:)\d*[1-9]+\d*\s*:\s*(?=\d{1,18}\s*,)\d*[1-9]+\d*\s*,\s*(t|T)\s*(?=\d{1,18}\s*:)\d*[1-9]+\d*\s*:\s*(?=\d{1,18}\s*$)\d*[1-9]+\d*\s*$)"))) {
                str.append_front("The units_format string passed, \"").append_back("\", does not match the format "
                                                                                   "required:\n\"Ma:b,Dc:x,Ty:z\", "
                                                                                   "where 'a', 'b', 'c', 'x', 'y', and "
                                                                                   "'z' represent any positive integer "
                                                                                   "number of 1-18 characters starting "
                                                                                   "with at least one non-zero "
                                                                                   "character.\n");
                throw std::invalid_argument(str.c_str());
            }
            str.strip("MmDdTt ");
            size_t colon_index = str.find(':');
            size_t comma_index = str.find(',');
            uint64_t m_denom = to_ull(str.substr(0, colon_index)); // denominator
            uint64_t m_num = to_ull(str.substr(colon_index + 1, comma_index)); // numerator
            str.erase_chars(0, comma_index + 1);
            colon_index = str.find(':');
            comma_index = str.find(',');
            uint64_t d_num = to_ull(str.substr(0, colon_index));
            uint64_t d_denom = to_ull(str.substr(colon_index + 1, comma_index));
            str.erase_chars(0, comma_index + 1);
            colon_index = str.find(':');
            comma_index = str.find(',');
            uint64_t t_denom = to_ull(str.substr(0, colon_index));
            uint64_t t_num = to_ull(str.substr(colon_index + 1, comma_index));
            G = 66'743; // 10^(-15) m^3 kg^-1 s^-2
            G *= m_num*d_num*d_num*d_num*t_num*t_num;
            G /= m_denom*d_denom*d_denom*d_denom*t_denom*t_denom;
            G /= 1'000'000'000'000'000; // correcting for the original G being 10^(15) times larger than it should be
            *G_vals++ = m_num;
            *G_vals++ = m_denom;
            *G_vals++ = d_num;
            *G_vals++ = d_denom;
            *G_vals++ = t_num;
            *G_vals = t_denom;
            G_vals -= 5;
        }
        void set_set() requires (collisions != 0) {
            if constexpr (memFreq)
                del_bods_m = new std::set<std::pair<ull_t, bod_t>, decltype(pair_bods_func)>{pair_bods_func};
            else
                del_bods_n = new std::set<bod_t, std::less<>>;
        }
        void clear_bodies() requires (memFreq != 0) { // makes sure the trajectories of all bodies are deleted
            // for (bod_t &bod : bods)
            //     bod.clear();
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.clear();});
        }
        void clear_bodies(vec_size_t &as_of) requires (memFreq != 0) {
            // if (as_of >= bods.size())
            //     return;
            // vec_size_t size = bods.size();
            // while (as_of < size)
            //     bods[as_of++].clear();
            std::for_each(bods.begin() + as_of, bods.end(), [this](bod_t &bod){bod.clear();});
        }
        void clear_bodies(vec_size_t &&as_of) requires (memFreq != 0) {
            clear_bodies(as_of);
        }
        void check_overlap() {
            bod_t *outer_b;
            bod_t *inner_b;
            num_bods = bods.size();
            for (outer = 0; outer < num_bods; ++outer) {
                outer_b = bods.data() + outer;
                inner_b = outer_b + 1;
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    if constexpr (mergeOverlappingBodies)
                        if (inner == outer) { // will be true at some point if there has been a merger
                            ++inner_b;
                            continue;
                        }
                    if ((outer_b->curr_pos - inner_b->curr_pos).magnitude() < outer_b->radius + inner_b->radius) {
                        if constexpr (!mergeOverlappingBodies) {
                            String str = "The bodies with id=";
                            str.append_back(outer_b->id).append_back(" and id=").append_back(inner_b->id);
                            str.append_back(" that were added to this system object overlap.\n");
                            throw overlapping_bodies_error(str.c_str());
                        }
                        *outer_b += *inner_b; // merges the two overlapping bodies
                        bods.erase(bods.begin() + inner); // thus the number of total bodies is reduced by 1
                        if (inner < outer) {
                            --outer;
                            --outer_b;
                        }
                        inner = 0; // have to recalculate possible mergers for newly created body
                        inner_b = bods.data();
                        --num_bods;
                        continue;
                    }
                    ++inner_b;
                }
            }
        }
        void cumulative_acc_and_pe(bod_t &b1, bod_t &b2) requires (memFreq != 0) {
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
        void cumulative_pe(bod_t &b1, bod_t &b2) requires (memFreq != 0) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/r12.magnitude();
            b1.pe += pot_energy;
            b2.pe += pot_energy;
        }
        template <bool mem, bool file>
        void take_euler_step() {
            for (bod_t &bod : bods) {
                bod.curr_pos += bod.curr_vel*dt;
                bod.curr_vel += bod.acc*dt;
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        template <bool mem, bool file>
        void take_modified_euler_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                      &predicted_vals) {
            outer = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += half_dt*(bod.curr_vel + std::get<1>(predicted_vals[outer]));
                bod.curr_vel += half_dt*(bod.acc + std::get<2>(predicted_vals[outer++]));
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        template <bool mem, bool file>
        void take_midpoint_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                &predicted_vals) {
            outer = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += dt*std::get<1>(predicted_vals[outer]);
                bod.curr_vel += dt*std::get<2>(predicted_vals[outer++]);
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        void create_energy_vectors() requires (memFreq != 0) {
            uint64_t iters_p1 = (iterations + 1)/memFreq;
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
            MEM_LOOP // zeros out total_pe, total_ke and pe of each body if memFreq, and always zeros all accelerations
            num_bods = bods.size();
            if constexpr (!memFreq) {
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner)
                        this->cumulative_acc(ref, bods[inner]);
                }
            }
            else {
                if (steps % memFreq) {
                    for (outer = 0; outer < num_bods; ++outer) {
                        bod_t &ref = bods[outer];
                        for (inner = outer + 1; inner < num_bods; ++inner)
                            this->cumulative_acc(ref, bods[inner]);
                    }
                    return;
                }
                else {
                    for (outer = 0; outer < num_bods; ++outer) {
                        bod_t &ref = bods[outer];
                        for (inner = outer + 1; inner < num_bods; ++inner)
                            this->cumulative_acc_and_pe(ref, bods[inner]);
                        ref.pe /= T{2};
                        pe[ref.id].push_back(ref.pe);
                        energy[ref.id].push_back(ref.curr_ke + ref.pe);
                        total_pe += ref.pe;
                        total_ke += ref.curr_ke;
                    }
                }
            }
            if constexpr (memFreq) {
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        void leapfrog_kdk_acc_e_and_step() {
            MEM_LOOP
            static auto loop = [this]<bool mem, bool file, bool only_e = false>{
                num_bods = bods.size();
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner) {
                        if constexpr (mem)
                            this->cumulative_acc_and_pe(ref, bods[inner]);
                        else if constexpr (only_e)
                            this->cumulative_pe(ref, bods[inner]);
                        else
                            this->cumulative_acc(ref, bods[inner]);
                    }
                    if constexpr (only_e) {
                        ref.pe /= T{2};
                        pe[ref.id].push_back(ref.pe);
                        energy[ref.id].push_back(ref.curr_ke + ref.pe);
                        total_pe += ref.pe;
                        total_ke += ref.curr_ke;
                    }
                    else {
                        /* KICK for half a step */
                        ref.curr_vel += half_dt*ref.acc;
                        if constexpr (mem) {
                            ref.pe /= T{2};
                            ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
                            pe[ref.id].push_back(ref.pe);
                            energy[ref.id].push_back(ref.curr_ke + ref.pe);
                            total_pe += ref.pe;
                            total_ke += ref.curr_ke;
                        }
                        if constexpr (file) {
                            // WRITE TO FILE
                        }
                    }
                }
            };
            FUNC_TEMPL_SELECT(loop.template operator(), return;)
            /* repeatedly evaluating the same if constexpr conditions DOES NOT MATTER (apart from slightly increasing
             * compilation time) given that they are evaluated at compile time, so there is never any runtime overhead*/
            if constexpr (collisions == overlap_coll_check) {
                // std::cout << "never get here!" << std::endl;
                loop.template operator()<false, false, true>();
            }
            // for (outer = 0; outer < num_bods; ++outer) {
            //     bod_t &ref = bods[outer];
            //     for (inner = outer + 1; inner < num_bods; ++inner)
            //         cumulative_acc_and_pe(ref, bods[inner]);
            //     ref.pe /= T{2};
            //     /* KICK for half a step */
            //     ref.curr_vel += half_dt*ref.acc;
            //     ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
            //     pe[ref.id].push_back(ref.pe);
            //     energy[ref.id].push_back(ref.curr_ke + ref.pe);
            //     total_pe += ref.pe;
            //     total_ke += ref.curr_ke;
            // }
            if constexpr (memFreq) { // in the FUNC_TEMPL_SELECT macro it is checked whether not steps % memFreq
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        void leapfrog_dkd_acc_e_and_step() {
            MEM_LOOP
            num_bods = bods.size();
            static auto loop = [this]<bool mem, bool file>{
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner)
                        cumulative_acc(ref, bods[inner]);
                    /* KICK for a full step */
                    ref.curr_vel += dt*ref.acc;
                    /* DRIFT for half a step */
                    ref.curr_pos += half_dt*ref.curr_vel;
                    /* the reason it is possible to update the outer body's position within the outer loop (without
                     * affecting the synchronisation of the particles) is because, by here, its effect (at its
                     * now-previous position) on all the other particles in the system has been calculated (all
                     * subsequent updates of the accelerations of the other particles no longer depend on the position
                     * of the outer body) */
                    if (mem) {
                        ref.add_pos_vel_ke(); // store new particle position, velocity and kinetic energy
                        total_ke += ref.curr_ke;
                    }
                    if (file) {
                        // WRITE TO FILE
                    }
                }
            };
            FUNC_TEMPL_SELECT(loop.template operator(), return;)
            /* a second loop is required to compute the potential energies based on the updated positions */
            if constexpr (memFreq) { // again, only reached if not steps % memFreq (checked within FUNC_TEMPL_SELECT)
                for (outer = 0; outer < num_bods; ++outer) {
                    bod_t &ref = bods[outer];
                    for (inner = outer + 1; inner < num_bods; ++inner)
                        cumulative_pe(ref, bods[inner]);
                    ref.pe /= T{2};
                    pe[ref.id].push_back(ref.pe);
                    energy[ref.id].push_back(ref.curr_ke + ref.pe);
                    total_pe += ref.pe;
                    if constexpr (collisions == overlap_coll_check)
                        total_ke += ref.curr_ke;
                }
                tot_pe.push_back(total_pe);
                tot_ke.push_back(total_ke);
                tot_e.push_back(total_pe + total_ke);
            }
        }
        void calc_energy() requires (memFreq != 0) {
            num_bods = bods.size();
            for (bod_t &bod : bods)
                bod.pe = T{};
            total_pe = total_ke = T{};
            for (outer = 0; outer < num_bods; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < num_bods; ++inner) {
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
        template <bool mem, bool file>
        void s_coll() { // "simple" (ahem, ahem) collision detection and evolution
            static std::map<long double, std::tuple<bod_t*, decltype(M{}*T{}), decltype(M{}*T{}), vector3D<T>>>
                    overlapping;
            outer = 0;
            num_bods = bods.size();
            bod_t *bod_o;
            bod_t *bod_i;
            bod_t *merging_bod;
            R rad_dist{};
            mom_t max_mom{min_tot_com_mom};
            mom_t curr_mom{};
            // long double min_dist = HUGE_VALL; // usually expands to infinity
            mom_t axis_mom_o;
            mom_t axis_mom_i;
            vec_size_t merging_index;
            bool first_merged;
            bod_t *merged;
            alignas(bod_t) char fake_body[sizeof(bod_t)];
            static std::set<uint64_t> merged_ids;
            merged_ids.clear();
            while (outer < num_bods) {
                first_merged = false;
                merging_bod = nullptr;
                overlapping.clear();
                bod_o = bods.data() + outer;
                bod_i = bod_o + 1; // to point to body just after bod_o
                for (inner = outer + 1; inner < num_bods; ++inner, ++bod_i) {
                    // bod_i = bods.data() + inner;
                    inner_loop:
                    if (inner == outer)
                        continue;
                    rad_dist = bod_o->radius + bod_i->radius;
                    vector3D<T> &&r12 = bod_i->curr_pos - bod_o->curr_pos;
                    long double &&dist = r12.magnitude(); // DEAL WITH ZERO DISTANCE CASE
                    r12.x /= dist; r12.y /= dist; r12.z /= dist; // more efficient than calling normalise()
                    if (rad_dist > dist) {
                        vector3D<T> &&com_vel = vel_com(bod_o, bod_i);
                        axis_mom_o = bod_o->momentum()*r12;
                        axis_mom_i = bod_i->momentum()*r12;
                        if((curr_mom = axis_mom_o - axis_mom_i - (bod_o->mass_ - bod_i->mass_)*(com_vel*r12))>=max_mom){
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
                    // bod_t &&merged = *bod_o + *merging_bod; // create body that is the result of the merger
                    if (!first_merged) {
                        merged = new(fake_body) bod_t{bod_o, merging_bod};//create body that is the result of the merger
                        if constexpr (memFreq) {
                            del_bods_m->emplace(steps, std::move(*bod_o));
                        }
                        else {
                            del_bods_n->emplace(std::move(*bod_o));
                        }
                        // bods.erase(bods.begin() + outer); // erase the merged bodies from the std::vector
                        // bod_o->~body();
                        bod_o = merged;
                        first_merged = true;
                    } else {
                        *bod_o += *merging_bod;
                    }
                    // must check that the merging body is not one from a previous merger at this time-step:
                    if (!merged_ids.contains(merging_bod->id)) {
                        if constexpr (memFreq) {
                            del_bods_m->emplace(steps, std::move(*merging_bod));
                        }
                        else {
                            del_bods_n->emplace(std::move(*merging_bod));
                        }
                    }
                    bods.erase(bods.begin() + merging_index);// - 1); // -1 because outer body was just deleted
                    // new position of primary merged body will be one less if merging body comes before
                    if (merging_index < outer) {
                        --outer;
                    }
                    --num_bods; // there is a net loss of 1 body
                    inner = 0;
                    bod_i = bods.data();
                    overlapping.clear();
                    merging_bod = nullptr;
                    goto inner_loop;
                }
                if (overlapping.size()) {
                    T o_vel;
                    T i_vel;
                    T o_minus_i;
                    T new_o_vel;
                    T new_i_vel;
                    long double avg_rest;
                    for (const auto &[_, tup] : overlapping) {
                        o_vel = std::get<1>(tup)/bod_o->mass_; // recalculating is cheaper than adding them to the tuple
                        i_vel = std::get<2>(tup)/bod_i->mass_; // up above
                        o_minus_i = o_vel - i_vel; // v1 - v2
                        avg_rest = (bod_o->rest_c + std::get<0>(tup)->rest_c)/2;
                        if (o_vel > 0 || i_vel < 0) {
                            if (o_vel <= 0 && i_vel < 0)
                                if (o_vel < i_vel) // bodies are already separating
                                    continue; // case for bodies having passed through each other
                            if (o_vel > 0 && i_vel >= 0)
                                if (o_vel < i_vel)
                                    continue; // case for bodies having passed through each other
                            new_o_vel = (std::get<1>(tup) + std::get<2>(tup) -
                                         std::get<0>(tup)->mass_*avg_rest*o_minus_i)/(bod_o->mass_ +
                                                                                      std::get<0>(tup)->mass_);
                            new_i_vel = (std::get<1>(tup) + std::get<2>(tup) +
                                         bod_o->mass_*avg_rest*o_minus_i)/(bod_o->mass_ + std::get<0>(tup)->mass_);
                            bod_o->curr_vel += (new_o_vel - o_vel)*std::get<3>(tup);
                            std::get<0>(tup)->curr_vel += (new_i_vel - i_vel)*std::get<3>(tup);
                        }
                    }
                }
                if (first_merged) {
                    // bod_o->rec_counter = steps;
                    // move the merged bodies into the deleted bodies set:
                    // bods.emplace_back(std::move(merged)); // move the new body into the std::vector storing all bodies
                    merged_ids.insert(bod_o->id);
                    bods[outer].~body(); // have to manually call dtor of body being moved into
                    copy(bods.data() + outer, fake_body, sizeof(bod_t)); // bypass constructors
                    // bod_o->~body_counter();
                    // if constexpr (file)
                    //     bod_o = bods.data() + outer;
                }
                // if constexpr (file) {
                //     // WRITE TO FILE
                // }
                ++outer;
            }
            if constexpr (memFreq || file) { // repeated constexpr conditionals don't make it to runtime!!
                // the repeated loop below is to avoid re-checking the conditions repeatedly inside the loop
                if (!merged_ids.empty()) {
                    // printf("Merged IDs:\n");
                    // for (const auto &id : merged_ids)
                    //     std::cout << id << std::endl;
                    for (bod_t &b : bods) {
                        if constexpr (mem) {
                            b.add_pos_vel_ke();
                        }
                        if constexpr (file) {
                            // WRITE TO FILE
                        }
                    }
                    return;
                }
                if (merged_ids.size() < pe.size()) {
                    for (bod_t &b : bods) {
                        if constexpr (mem) {
                            b.add_pos_vel_ke();
                        }
                        // bod_o->rec_counter = steps;
                        if constexpr (memFreq) {
                            if (!merged_ids.contains(b.id)) {
                                pe.emplace(b.id, std::vector<T>{}); // add a std::vector to store its potential energies
                                energy.emplace(b.id, std::vector<T>{}); // same for its total energy (stores its own KE)
                            }
                        }
                        if constexpr (file) {
                            // WRITE TO FILE
                        }
                    }
                    return;
                }
                for (bod_t &b : bods) {
                    if constexpr (mem) {
                        b.add_pos_vel_ke();
                    }
                    // bod_o->rec_counter = steps;
                    if constexpr (memFreq) {
                        pe.emplace(b.id, std::vector<T>{}); // add a std::vector to store its potential energies
                        energy.emplace(b.id, std::vector<T>{}); // same for its total energy (it stores its own KE)
                    }
                    if constexpr (file) {
                        // WRITE TO FILE
                    }
                }
            }
        }
    public:
        void analysis() requires (memFreq != 0) { // debugging method - mark for removal
            printf("---------------ANALYSIS--------------\nBods:\n");
            for (const auto& bod : bods) {
                std::cout << "\nBody: " << bod <<
                          "\nPositions size: " << bod.positions.size() <<
                          ", Velocity size: " << bod.velocities.size() <<
                          "\nKE size: " << bod.energies.size() << ", PE size: " << pe[bod.id].size() <<
                          ", Energy size: " << energy[bod.id].size() << '\n' << std::endl;
            }
            printf("-----------\nDel_Bods:\n");
            for (const auto& [it, bod] : *del_bods_m) {
                std::cout << "Iteration: " << it << "\nBody: " << bod <<
                "\nPositions size: " << bod.positions.size() <<
                ", Velocity size: " << bod.velocities.size() <<
                "\nKE size: " << bod.energies.size() << ", PE size: " << pe[bod.id].size() <<
                ", Energy size: " << energy[bod.id].size() << '\n' << std::endl;
            }
        }
        std::vector<bod_t>& get_bods() {
            return this->bods;
        }
        std::set<std::pair<uint64_t, bod_t>, decltype(pair_bods_func)>& get_del_bods() requires (memFreq != 0) {
            return *this->del_bods_m;
        }
    private:
        static inline uint64_t nsys_file_size(uint64_t n_bods) noexcept {
            char rem = (char) ((n_bods*sizeof(nsys_chunk)) % 4);
            return sizeof(nsys_header) + n_bods*sizeof(nsys_chunk) + (rem ? 4 - rem : 0);
        }
        // path is guaranteed not to be nullptr here:
        void load_from_nsys(const char *path, bool alloc_chunks) requires (std::is_fundamental_v<M> &&
                                                                           std::is_fundamental_v<R> &&
                                                                           std::is_fundamental_v<T> && CHAR_BIT == 8) {
#ifndef _WIN32
            /* Unfortunately, there is absolutely no POSIX-conforming standard way of ensuring that files larger than
             * 2GB can be worked with. The size of off_t can be checked to ensure it is 64 bits wide, but nothing can be
             * done if it is not. There are many non-standard extensions, such as defining the _FILE_OFFSET_BITS macro
             * as 64 (as a GNU extension), but as these are not part of the standard, I do not include them. */
            struct stat buffer{};
            if (stat(path, &buffer) == -1)
                throw nsys_load_error{"Error obtaining .nsys file information. No data loaded.\n"};
            if (!S_ISREG(buffer.st_mode))
                throw nsys_load_error{".nsys file provided is not a regular file. Could not load data.\n"};
            if (buffer.st_size < sizeof(nsys_header))
                throw nsys_load_error{"Insufficient .nsys file size. No data loaded.\n"};
#else
            WIN32_FILE_ATTRIBUTE_DATA buffer{};
            if (!GetFileAttributesExA(path, GetFileExInfoStandard, &buffer))
                throw nsys_load_error{"Error obtaining .nsys file attributes. No n-body data loaded.\n"};
            if (buffer.dwFileAttributes != FILE_ATTRIBUTE_NORMAL)
                throw nsys_load_error{".nsys file provided is not a regular file. No data loaded.\n"};
            uint64_t fileSize = buffer.dwFileSizeLow + (buffer.dwFileSizeHigh << 32);
            if (fileSize < sizeof(nsys_header))
                throw nsys_load_error{".nsys file size insufficient. Could not load data."};
#endif
            istream = new std::ifstream{path, std::ios_base::in | std::ios_base::binary};
            if (!*istream) {
                delete istream;
                throw nsys_load_error{"Error opening .nsys file. Could not load data.\n"};
            }
            nsys_header header;
            istream->read((char *) &header, sizeof(nsys_header));
            if (header.signature[0] != 'N' && header.signature[1] != 'S') {
                istream->close();
                delete istream;
                throw nsys_load_error{"Invalid .nsys format: invalid file signature.\n"};
            }
            if ((header.d_types & 0b1000000) != 0) { // 0b10000000 == 128, but I use the binary repr. for clarity
                istream->close();
                delete istream;
                throw nsys_load_error{"Invalid .nsys format: msb occupied in d_types field (should be zero).\n"};
            } // the following macro saves about 60 lines of code
#define NSYS_TYPE_CHECK(type, bin_num1, bin_num2, str) \
            if constexpr (std::floating_point<type>) { \
                if (!(header.d_types & bin_num1)) {    \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " floating-point format, but template parameter "#type" is integral.\n"}; \
                } \
            } \
            else if constexpr (std::signed_integral<type>) { \
                if (header.d_types & bin_num1) {       \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " floating-point format, but template parameter "#type" is integral.\n"}; \
                }/* I have a separate conditional for the case below so that a different error message can be logged */\
                if (!(header.d_types & bin_num2)) {    \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " unsigned integral format, but template parameter "#type" is signed " \
                                          "integral.\n"}; \
                } \
            } \
            else { \
                if (header.d_types & bin_num1) {       \
                    istream->close();\
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          " floating-point format, but template parameter "#type" is integral.\n"}; \
                } \
                if (header.d_types & bin_num2) {       \
                    istream->close(); \
                    delete istream; \
                    throw nsys_load_error{"Invalid .nsys format: d_types field states that "#str" values are stored in"\
                                          "signed integral format, but template parameter "#type" is unsigned " \
                                          "integral.\n"}; \
                } \
            } \
            if (header.type##_size != sizeof(type)) {  \
                istream->close(); \
                delete istream; \
                char error_msg[128]; \
                snprintf(error_msg, 128, "Size mismatch error: size of "#type" template parameter (%d bytes) does " \
                                         "not match its reported size (%d bytes) in the .nsys file.\n", \
                        sizeof(type), header.type##_size); /* includes null terminator */ \
                throw nsys_load_error{error_msg}; \
            }
            NSYS_TYPE_CHECK(M, 0b00000010, 0b00010000, mass)
            NSYS_TYPE_CHECK(R, 0b00000100, 0b00100000, radius)
            NSYS_TYPE_CHECK(T, 0b00001000, 0b01000000, position and velocity)
#undef NSYS_TYPE_CHECK
            if (header.rest_coeff_size > sizeof(long double)) {
                istream->close();
                delete istream;
                char error_msg[156];
                snprintf(error_msg, 156, "Size mismatch error: coefficient of restitution data size (%d bytes) larger "
                                         "than widest floating point data type available (long double = %d bytes).\n",
                         header.rest_coeff_size, sizeof(long double)); // includes null terminator
                throw nsys_load_error{error_msg};
            }
            void (*rest_c_func)(long double *, const char *) = +[](long double *coeff, const char *ptr) {
                char *coeff_ptr = (char *) coeff;
                for (unsigned char counter = 0; counter < sizeof(long double); ++counter)
                    *coeff_ptr++ = *ptr++;
            };
            if (header.rest_coeff_size != sizeof(long double)) {
                if (header.rest_coeff_size == sizeof(double)) {
                    rest_c_func = +[](long double *coeff, const char *ptr) {
                        static double coeff_d;
                        char *coeff_ptr = (char *) &coeff_d;
                        for (unsigned char counter = 0; counter < sizeof(double); ++counter)
                            *coeff_ptr++ = *ptr++;
                        *coeff = coeff_d;
                    };
                } else if (header.rest_coeff_size == sizeof(float)) {
                    rest_c_func = +[](long double *coeff, const char *ptr) {
                        static float coeff_f;
                        char *coeff_ptr = (char *) &coeff_f;
                        for (unsigned char counter = 0; counter < sizeof(float); ++counter)
                            *coeff_ptr++ = *ptr++;
                        *coeff = coeff_f;
                    };
                } else {
                    istream->close();
                    delete istream;
                    char error_msg[192];
                    snprintf(error_msg, 192, "Size mismatch error: coefficient of restitution data size (%d bytes) "
                                             "does not match any floating point data type available:\n"
                                             "float = %zu bytes,\n"
                                             "double = %zu bytes,\n"
                                             "long double = %zu bytes.\n",
                             header.rest_coeff_size, sizeof(float), sizeof(double), sizeof(long double));
                    throw nsys_load_error{error_msg};
                }
            }
#ifndef _WIN32
            if (buffer.st_size != nsys_file_size(header.num_bodies)) {
#else
            if (fileSize != nsys_file_size(header.num_bodies)) {
#endif
                istream->close();
                delete istream;
                char error_msg[196];
                snprintf(error_msg, 196, "Error: file size does not match expected file size based on number of bodies."
                                         "\nExpected: %llu bytes. Instead got: %zu bytes. No data loaded.\n",
#ifndef _WIN32
                         (unsigned long long) nsys_file_size(header.num_bodies), (size_t) buffer.st_size);
#else
                         (unsigned long long) nsys_file_size(header.num_bodies), (size_t) fileSize);
#endif
                throw nsys_load_error{error_msg};
            }
            /* For errors with the size of the data type of the time value, I do not throw any exception, as the time
             * value is not needed for the functioning of the gtd::system<> class. Nonetheless, if there were errors
             * present in its format, time_elapsed is set to the largest possible value of a long double. This can be
             * queried via the gtd::system<>::elapsed_time member function. */
            if (header.d_types & 0b00000001) { // case for floating-point time value
                switch (header.time_size) {
                    case sizeof(long double):
                        copy(&this->time_elapsed, header.time_point, sizeof(long double));
                        copy(&this->dt, header.delta_t, sizeof(long double));
                        break;
                    case sizeof(double): {
                        double time_val;
                        copy(&time_val, header.time_point, sizeof(double));
                        this->time_elapsed = time_val;
                        copy(&time_val, header.delta_t, sizeof(double));
                        this->dt = time_val;
                        break;
                    }
                    case sizeof(float): { // almost certainly == 4
                        float time_val;
                        copy(&time_val, header.time_point, sizeof(float));
                        this->time_elapsed = time_val;
                        copy(&time_val, header.delta_t, sizeof(float));
                        this->dt = time_val;
                        break;
                    }
                    default:
                        istream->close();
                        delete istream;
                        char error_msg[208];
                        snprintf(error_msg, 208, "Size mismatch error: reported size of data type for time values in "
                                                 ".nsys file (%d bytes) larger than any floating point data type "
                                                 "available:\n"
                                                 "float = %zu bytes,\n"
                                                 "double = %zu bytes,\n"
                                                 "long double = %zu bytes.\n",
                                 header.time_size, sizeof(float), sizeof(double), sizeof(long double));
                        throw nsys_load_error{error_msg};
                        // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                        //                      std::numeric_limits<long double>::quiet_NaN() :
                        //                      std::numeric_limits<long double>::max(); // error case
                }
            }
            else { // case for integral time value
                if constexpr (sizeof(unsigned long long) >= sizeof(uint64_t)) {
                    if (header.time_size > sizeof(unsigned long long)) {
                        // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                        //                      std::numeric_limits<long double>::quiet_NaN() :
                        //                      std::numeric_limits<long double>::max(); // error case
                    }
                    unsigned long long time_val = 0;
                    copy(&time_val, header.time_point, header.time_size);
                    this->time_elapsed = time_val;
                    copy(&time_val, header.delta_t, header.time_size);
                    this->dt = time_val;
                }
                else {
                    if (header.time_size > sizeof(uint64_t)) {
                        // this->time_elapsed = std::numeric_limits<long double>::has_quiet_NaN ?
                        //                      std::numeric_limits<long double>::quiet_NaN() :
                        //                      std::numeric_limits<long double>::max(); // error case
                    }
                    uint64_t time_val = 0;
                    copy(&time_val, header.time_point, header.time_size);
                    this->time_elapsed = time_val;
                    copy(&time_val, header.delta_t, header.time_size);
                    this->dt = time_val;
                }
            }
            this->half_dt = this->dt/2;
            this->G = 66'743;
            this->G *= header.m_num*header.d_num*header.d_num*header.d_num*header.t_num*header.t_num;
            this->G /= header.m_denom*header.d_denom*header.d_denom*header.d_denom*header.t_denom*header.t_denom;
            this->G /= 1'000'000'000'000'000;
            *G_vals++ = header.m_num;
            *G_vals++ = header.m_denom;
            *G_vals++ = header.d_num;
            *G_vals++ = header.d_denom;
            *G_vals++ = header.t_num;
            *G_vals = header.t_denom;
            G_vals -= 5;
            this->bods.clear();
            if (!header.num_bodies) {
                istream->close();
                delete istream;
                return;
            }
            this->bods.reserve(header.num_bodies); // avoids memory reallocation during the loop
            long double *restitution;
            if (alloc_chunks) {
                /* Option 1: read in all data from chunks into allocated memory, and construct bodies in the "bods"
                 * std::vector using the data from the chunks. Much faster than option 2, but uses much more memory. */
                nsys_chunk *chunks = new nsys_chunk[header.num_bodies];
                nsys_chunk *chunk = chunks;
                istream->read((char *) chunks, header.num_bodies*sizeof(nsys_chunk));
                istream->close();
                delete istream;
                while (header.num_bodies --> 0) {
                    rest_c_func(&restitution, chunk->r_coeff);
                    /* I create a new, nasty constructor in the gtd::body<> class which accepts individual components of
                     * the position and velocity (rather than two vectors) to avoid making any unnecessary copies. */
                    this->bods.emplace_back(chunk->mass, chunk->radius, chunk->xpos, chunk->ypos, chunk->zpos,
                                            chunk->xvel, chunk->yvel, chunk->zvel, restitution);
                    ++chunk;
                }
                delete [] chunks;
                return;
            }
            /* Option 2: read data in from the file chunk by chunk, constructing each gtd::body<> within the "bods"
             * std::vector after each read. Slower than option 1 as requires N times more I/O function calls (where N is
             * number of bodies), but reduces memory footprint. */
            nsys_chunk chunk;
            while (header.num_bodies --> 0) {
                istream->read((char *) &chunk, sizeof(nsys_chunk));
                rest_c_func(&restitution, chunk.r_coeff);
                this->bods.emplace_back(chunk.mass, chunk.radius, chunk.xpos, chunk.ypos, chunk.zpos,
                                        chunk.xvel, chunk.yvel, chunk.zvel, restitution);
            }
            istream->close();
            delete istream;
        }
        void check_dt() const {
            if (this->dt < 0)
                throw std::invalid_argument{"Time-step cannot be negative.\n"};
        }
        static inline bool check_option(int option) noexcept {
            uint16_t loword = option & 0x0000ffff;
            uint16_t hiword = option >> 16;
            return !(loword & (loword - 1)) || !loword || loword > rk4 ||
                   !(hiword & (hiword - 1)) || hiword > 2;
        }
        // constexpr void check_coll() {
        //     static_assert(collisions <= pred_coll_check && !(collisions & (collisions + 1)),
        //                   "Invalid collision-checking option.");
        // }
        void print_progress() const noexcept requires (prog) {
#ifndef _WIN32
            printf(CYAN_TXT_START "Iteration " BLUE_TXT_START "%llu" RED_TXT_START "/" MAGENTA_TXT_START
                   "%llu\r", this->steps, this->iterations);
#else
            printf("Iteration %llu/%llu\r", step, iterations);
#endif
        }
        void print_conclusion(const std::chrono::time_point<std::chrono::high_resolution_clock> &start,
                              std::chrono::nanoseconds &total) requires (prog) {
            total = std::chrono::high_resolution_clock::now() - start;
            std::cout << RESET_TXT_FLAGS << BLACK_TXT("\n--------------------Done--------------------\n") <<
            UNDERLINED_TXT_START GREEN_TXT_START << this->method_str() <<
            RESET_TXT_FLAGS WHITE_TXT(" - time elapsed: ") BOLD_TXT_START YELLOW_TXT_START << total.count()/BILLION <<
            RESET_TXT_FLAGS WHITE_TXT_START " second" << "s"[total == std::chrono::seconds{1}] <<
            RESET_TXT_FLAGS << std::endl;
        }
    public:
        static_assert(collisions <= pred_coll_check && !(collisions & (collisions + 1)),
                      "Invalid collision-checking option.");
        typedef struct nsys_file_header {
            const char signature[2] = {'N', 'S'};
            const char d_types = 1 + (std::floating_point<M> << 1) +
                                     (std::floating_point<R> << 2) +
                                     (std::floating_point<T> << 3) +
                                     (std::signed_integral<M> << 4) +
                                     (std::signed_integral<R> << 5) +
                                     (std::signed_integral<T> << 6);
            const char time_size = sizeof(long double);
            char time_point[16]{};
            char delta_t[16]{};
            const unsigned char M_size = sizeof(M);
            const unsigned char R_size = sizeof(R);
            const unsigned char T_size = sizeof(T);
            const unsigned char rest_coeff_size = sizeof(long double);
            uint64_t m_num{};
            uint64_t m_denom{};
            uint64_t d_num{};
            uint64_t d_denom{};
            uint64_t t_num{};
            uint64_t t_denom{};
            uint64_t num_bodies{};
        } nsys_header; // size == 96 bytes
        typedef struct nsys_file_chunk {
            char r_coeff[16];
            M mass{};
            R radius{};
            T xpos{};
            T ypos{};
            T zpos{};
            T xvel{};
            T yvel{};
            T zvel{};
        } nsys_chunk;
        system() : G{G_SI} {this->check_coll(); if constexpr (collisions) this->set_set();}
        explicit system(const char *nsys_path, bool allocate_chunks = false) {
            if (nsys_path == nullptr)
                throw std::invalid_argument{"nullptr passed as .nsys file path.\n"};
            this->load_from_nsys(nsys_path);
            this->check_overlap();
            if constexpr (collisions)
                this->set_set();
        }
        system(const std::initializer_list<bod_t> &list) : bods{list}, G{G_SI} {
            check_overlap();
            // check_coll();
            if constexpr (memFreq)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        system(std::initializer_list<bod_t> &&list) : bods{std::move(list)}, G{G_SI} {
            check_overlap();
            // check_coll();
            if constexpr (memFreq)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        explicit system(long double timestep, uint64_t num_iterations, const char *units_format) :
                dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            // check_coll();
            this->check_dt();
            parse_units_format(units_format);
            if constexpr (collisions)
                this->set_set();
        }
        template <uint64_t mF>
        system(const std::vector<body<M, R, T, mF>> &bodies, long double timestep = 1,
               uint64_t num_iterations = 1'000, const char *units_format = "M1:1,D1:1,T1:1") :
               bods{bodies.begin(), bodies.end()}, dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            this->check_dt();
            check_overlap();
            // check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq && mF) // no need to clear bodies if the ones passed do not record history
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        system(std::vector<bod_t> &&bodies, long double timestep = 1,
               uint64_t num_iterations = 1'000, const char *units_format = "M1:1,D1:1,T1:1") :
               bods{std::move(bodies)}, dt{timestep}, iterations{num_iterations} {
            this->check_dt();
            check_overlap();
            // check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        template <uint64_t mF>
        system(std::vector<body<M, R, T, mF>> &&bodies, long double timestep = 1,
               uint64_t num_iterations = 1'000, const char *units_format = "M1:1,D1:1,T1:1") :
               dt{timestep}, iterations{num_iterations} {
            this->check_dt();
            std::for_each(bodies.begin(), bodies.end(),
                          [this](body<M, R, T, mF> &b){bods.emplace_back(std::move(b));});
            check_overlap();
            // check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq && mF)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        /* copy constructors are made to only copy the variables seen below: it does not make sense for a new system
         * object (even when copy-constructed) to be "evolved" if it has not gone through the evolution itself */
        template <bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        // so a system can be constructed from another without same checks
        system(const system<M, R, T, prg, mrg, coll, mF, fF, bF> &other) :
               bods{other.bods.begin(), other.bods.end()}, dt{other.dt}, iterations{other.iterations}, G{other.G} {
            check_overlap();
            // check_coll();
            if constexpr (memFreq && mF)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        template <bool prg, bool mrg, int coll, uint64_t fF, bool bF>
        system(system<M, R, T, prg, mrg, coll, memFreq, fF, bF> &&other) :
               bods{std::move(other.bods)}, dt{other.dt}, iterations{other.iterations}, G{other.G} {
            check_overlap();
            // check_coll();
            if constexpr (memFreq)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        template <bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        system(system<M, R, T, prg, mrg, coll, mF, fF, bF> &&other) :
               dt{other.dt}, iterations{other.iterations}, G{other.G} {
            std::for_each(other.bods.begin(), other.bods.end(),
                          [this](body<M, R, T, mF> &b){bods.emplace_back(std::move(b));});
            check_overlap();
            // check_coll();
            if constexpr (memFreq && mF)
                clear_bodies();
            if constexpr (collisions)
                this->set_set();
        }
        vec_size_t num_bodies() {
            return bods.size();
        }
        long double elapsed_time() const noexcept {
            return this->time_elapsed; // will either be NaN or LDBL_MAX if bad value read in from .nsys file
        }
        bool set_iterations(uint64_t number) noexcept {
            if (!number) // cannot have zero iterations
                return false;
            iterations = number;
            return true;
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
        void set_G_units(const char *units_format) {
            this->parse_units_format(units_format);
        } // the size the .nsys file would be at the time the function is called:
        uint64_t nsys_file_size() requires (std::is_fundamental_v<M> &&
                                            std::is_fundamental_v<R> &&
                                            std::is_fundamental_v<T> && CHAR_BIT == 8) {
            return nsys_file_size(this->bods.size());
        }
        sys_t &add_body(const bod_t &bod) {
            bods.emplace_back(bod);
            if constexpr (memFreq)
                bods.back().clear(); // all bodies within a system object must start out without an evolution
            check_overlap();
            return *this;
        }
        sys_t &add_body(bod_t &&bod) {
            if constexpr (memFreq)
                bod.clear();
            bods.emplace_back(std::move(bod));
            check_overlap();
            return *this;
        }
        template <typename ...Args>
        sys_t &emplace_body(Args&& ...args) { // to allow a body to be constructed in-place, within the system object
            bods.emplace_back(std::forward<Args>(args)...);
            this->check_overlap();
            return *this;
        }
        sys_t &add_bodies(const std::vector<bod_t> &bodies) {
            if (!bodies.size())
                return *this;
            num_bods = bods.size();
            bods.insert(bods.end(), bodies.begin(), bodies.end());
            if constexpr (memFreq)
                clear_bodies(num_bods);
            check_overlap();
            return *this;
        }
        sys_t &add_bodies(std::vector<bod_t> &&bodies) {
            if (!bodies.size())
                return *this;
            for (auto &b : bodies) {
                if constexpr (memFreq)
                    b.clear();
                bods.emplace_back(std::move(b));
            }
            check_overlap();
            return *this;
        }
        const bod_t &back() const {
            if (this->bods.size())
                return this->bods.back();
            throw std::logic_error{"gtd::system<>::back cannot be called on an empty gtd::system object.\n"};
        }
        const bod_t &front() const {
            if (this->bods.size())
                return this->bods.front();
            throw std::logic_error{"gtd::system<>::front cannot be called on an empty gtd::system object.\n"};
        }
        const bod_t &get_body(uint64_t id) const {
            for (const auto &b : *this) // this is where it would be nice to be using a set!!
                if (b.id == id)
                    return b;
            throw std::invalid_argument("The id passed does not correspond to any body present in this system "
                                        "object.\n");
        }
        bool remove_body(uint64_t id) {
            auto end_it = bods.cend();
            for (typename std::vector<bod_t>::const_iterator it = bods.cbegin(); it < end_it; ++it)
                if ((*it).id == id) {
                    bods.erase(it);
                    return true;
                }
            return false;
        }
        void reset() requires (memFreq != 0) {
            // for (bod_t &bod : bods)
            //     bod.reset(true);
            if constexpr (collisions) {
                bods.insert(bods.end(), std::make_move_iterator(del_bods_m->begin()),
                            std::make_move_iterator(del_bods_m->end()));
                del_bods_m->clear();
            }
            std::for_each(bods.begin(), bods.end(), [this](bod_t &bod){bod.reset(true);});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
            time_elapsed = 0;
        }
        /* reset() deletes the evolution of the bodies AND returns them to their original positions, velocities, and
         * kinetic energies, whilst clear_evolution() deletes the evolution BUT retains the current positions,
         * velocities and kinetic energies of the bodies. */
        void clear_evolution() requires (memFreq != 0) {
            // for (bod_t &bod : bods)
            //     bod.clear();
            // bods.insert(bods.end(), std::make_move_iterator(del_bods_m->begin()),
            //             std::make_move_iterator(del_bods_m->end()));
            if constexpr (collisions)
                del_bods_m->clear();
            std::for_each(bods.begin(), bods.end(), [](bod_t &bod){bod.clear();});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
            time_elapsed = 0;
        }
    private:
        const char *&method_str() requires (prog) {
            static const char *ptr{};
            return ptr;
        }
    public:
        const std::chrono::nanoseconds &evolve(int integration_method = leapfrog_kdk) {
            static std::chrono::time_point<std::chrono::high_resolution_clock> start;
            static std::chrono::nanoseconds total;
            num_bods = bods.size();
            if (!num_bods || !check_option(integration_method))
                return total = std::chrono::nanoseconds::zero();
            // time_t start = time(nullptr);
            steps = 0;
            if constexpr (memFreq) {
                if (evolved)
                    this->clear_evolution();
                this->create_energy_vectors();
            }
#ifndef _WIN32
            if constexpr (prog)
                std::cout << BOLD_TXT_START;
#endif
            start = std::chrono::high_resolution_clock::now();
            if ((integration_method & two_body) == two_body) {
                if (num_bods != 2)
                    throw two_body_error();
                /* no force approximation techniques performed here, as the integration is just for two bodies */
            }
            else if ((integration_method & euler) == euler) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    while (steps < iterations) {
                        if constexpr (collisions == pred_coll_check) {

                        }
                        this->calc_acc_and_e();
                        ++steps;
                        // this->take_euler_step();
                        FUNC_TEMPL_SELECT(this->take_euler_step, EMPTY)
                        /* by having "prog" as a template parameter, it avoids having to re-evaluate whether it is true
                         * or not within the loop, since the "if constexpr ()" branches get rejected at compile-time */
                        // if constexpr (collisions == overlap_coll_check)
                        //     this->s_coll();
                        if constexpr (prog)
                            this->print_progress();
                    }
                    if constexpr (prog)
                        this->method_str() = "Euler";
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
                    while (steps < iterations) {
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
                        // this->take_modified_euler_step(predicted);
                        ++steps;
                        FUNC_TEMPL_SELECT(this->take_modified_euler_step, EMPTY, predicted)
                        for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup : predicted)
                            std::get<2>(tup).make_zero();
                        if constexpr (collisions == overlap_coll_check) {
                            // this->s_coll();
                            num_bods = bods.size(); // number of bodies can change through mergers
                        }
                        if constexpr (prog)
                            this->print_progress();
                    }
                    if constexpr (prog)
                        this->method_str() = "Modified Euler";
                }
            }
            else if ((integration_method & midpoint) == midpoint) {
                std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>> predicted{num_bods};
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    while (steps < iterations) {
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
                        // this->take_midpoint_step(predicted);
                        ++steps;
                        FUNC_TEMPL_SELECT(this->take_midpoint_step, EMPTY, predicted)
                        for (std::tuple<vector3D<T>, vector3D<T>, vector3D<T>> &tup : predicted)
                            std::get<2>(tup).make_zero();
                        if constexpr (collisions == overlap_coll_check) {
                            // this->s_coll();
                            num_bods = bods.size();
                        }
                        if constexpr (prog)
                            this->print_progress();
                    }
                    if constexpr (prog)
                        this->method_str() = "Midpoint method";
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
                        // if constexpr (collisions == overlap_coll_check)
                        //     this->s_coll();
                        if constexpr (prog)
                            this->print_progress();
                    }
                    if constexpr (prog)
                        this->method_str() = "Leapfrog KDK";
                }
                goto bye;
            }
            else if ((integration_method & leapfrog_dkd) == leapfrog_dkd) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    if constexpr (memFreq)
                        this->calc_energy();
                    while (steps++ < iterations) {
                        for (bod_t &bod : bods)
                            /* DRIFT for half a step */
                            bod.curr_pos += half_dt*bod.curr_vel;
                        /* update accelerations, then KICK for a full step and DRIFT for half a step, then update E */
                        this->leapfrog_dkd_acc_e_and_step();
                        // if constexpr (collisions == overlap_coll_check)
                        //     this->s_coll();
                        if constexpr (prog)
                            this->print_progress();
                    }
                    if constexpr (prog)
                        this->method_str() = "Leapfrog DKD";
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
            if constexpr (memFreq) {
                if (!(steps % memFreq)) {
                    this->calc_energy();
                }
            }
            bye:
            if constexpr (prog)
                print_conclusion(start, total);
            else
                total = std::chrono::high_resolution_clock::now() - start;
            evolved = true;
            prev_dt = dt;
            prev_iterations = iterations;
            time_elapsed = iterations*dt;
            return total;
        }
        std::pair<T*, long double> *pe_extrema() const requires (memFreq != 0) {
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
            min_max->second = min_index*prev_dt*memFreq; // assign times at which extrema occurred
            (min_max + 1)->second = max_index*prev_dt*memFreq;
            return min_max;
        }
        std::pair<T*, long double> *ke_extrema() const requires (memFreq != 0) {
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
            min_max->second = min_index*prev_dt*memFreq;
            (min_max + 1)->second = max_index*prev_dt*memFreq;
            return min_max;
        }
        std::pair<T*, long double> *tot_e_extrema() const requires (memFreq != 0) {
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
            min_max->second = min_index*prev_dt*memFreq;
            (min_max + 1)->second = max_index*prev_dt*memFreq;
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
        std::streamoff to_nsys(const char *path) const requires (std::is_fundamental_v<M> &&
                                                                 std::is_fundamental_v<R> &&
                                                                 std::is_fundamental_v<T> && CHAR_BIT == 8) {
            static nsys_header header;
            static nsys_chunk chunk;
            static std::ofstream::pos_type tell;
            if (path == nullptr)
                return 0;
            ostream = new std::ofstream{path, std::ios_base::trunc | std::ios_base::binary};
            if (!*ostream) {
                delete ostream;
                return 0;
            }
            copy(header.time_point, &this->time_elapsed, sizeof(long double));
            copy(header.delta_t, &this->dt, sizeof(long double));
            copy(&header.m_num, this->G_vals, 6*sizeof(uint64_t)); // I include sizeof just so intent is clear
            header.num_bodies = this->bods.size();
            ostream->write((char *) &header, sizeof(nsys_header));
            if (!header.num_bodies) { // case for no bodies - no reason to prohibit this
                ostream->close();
                delete ostream;
                return sizeof(nsys_header);
            }
            for (const bod_t &bod : bods) {
                // chunk.r_coeff = bod.rest_c;
                copy(chunk.r_coeff, &bod.rest_c, sizeof(long double));
                chunk.mass = bod.mass_;
                chunk.radius = bod.radius;
                chunk.xpos = bod.curr_pos.x;
                chunk.ypos = bod.curr_pos.y;
                chunk.zpos = bod.curr_pos.z;
                chunk.xvel = bod.curr_vel.x;
                chunk.yvel = bod.curr_vel.y;
                chunk.zvel = bod.curr_vel.z;
                ostream->write((char *) &chunk, sizeof(nsys_chunk));
            }
            tell = ostream->tellp();
            if (char rem = (char) (tell % 4); rem) {
                ostream->write("\0\0\0", 4 - rem); // align EOF with dword boundary - at most 3 bytes would be written
                tell += (std::streamoff) (4 - rem);
            }
            ostream->close();
            delete ostream;
            return tell;
        }
        // static std::pair<nsys_header, std::vector<bod_t>> load_nsys(const char *path, bool only_header = false) {
        //     if (path == nullptr)
        //         return {};
        // }
        bool write_trajectories(const String &path = def_path, bool binary = false) requires (memFreq != 0) {
            /* This method can only be called if memFreq is non-zero, or else there would be no history to output. */
            /* Writes the trajectories (historically) of all the bodies to a file, including all the energy values. For
             * text files, units are not written as these will depend entirely on the units_format string passed to the
             * system<M, R, T> object in its constructor. */
            if (!evolved || !bods.size())
                return false;
            if (binary) {
                if (&path == &def_path)
                    def_path.append_back(get_date_and_time()).append_back(".nbod");
                if (&path == &def_path)
                    def_path.erase_chars(def_path.get_length() - 30);
                return true;
            } // needs work - does NOT output energies
            if (&path == &def_path)
                def_path.append_back(get_date_and_time()).append_back(".csv");
            std::ofstream out(path.c_str(), std::ios_base::trunc);
            if (!out)
                return false;
            // unsigned long long count = 0;
            vec_size_t i;
            const ull_t num_max_els = (prev_iterations + (memFreq != 1))/memFreq; // CHECK THIS IS CORRECT
            // std::cout << "size of bods: " << bods.size() << ", size of del_bods: " << del_bods_m->size() << std::endl;
            if constexpr (collisions) {
                ull_t d_i;
                vec_size_t max_index;
                for (const auto &[destruction_it, bod] : *del_bods_m) {
                    max_index = bod.positions.size();
                    // d_i = (destruction_it + memFreq)/memFreq - max_index - 1;
                    if (destruction_it % memFreq) {
                        d_i = destruction_it/memFreq + 1 - max_index;
                    } else {
                        d_i = destruction_it/memFreq - max_index;
                    }
                    // d_i = destruction_it/memFreq;
                    // d_i += (max_index > d_i) - max_index;
                    // std::cout << "Positions size: " << bod.positions.size() << ", pe size: " << pe[bod.id].size()
                    // << std::endl;
                    std::cout << "body id: " << bod.id << std::endl << std::endl;
                    out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                           "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,"
                           "kinetic_energy,potential_energy,total_energy\r\n";
                    for (i = 0; i < max_index; ++i)
                        out << d_i++*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                            << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                            << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                            << energy[bod.id][i] << "\r\n";
                    out << ",,,,,,,,,\r\n";
                }
                for (const bod_t &bod : bods) {
                    max_index = bod.positions.size();
                    d_i = prev_iterations/memFreq - max_index + 1;
                    // std::cout << "max_index: " << max_index << ", d_i: " << d_i << std::endl;
                    // std::cout << "body id: " << bod.id << std::endl;
                    out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                           "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,"
                           "kinetic_energy,potential_energy,total_energy\r\n";
                    for (i = 0; i < max_index; ++i)
                        out << d_i++*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                            << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                            << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                            << energy[bod.id][i] << "\r\n";
                    out << ",,,,,,,,,\r\n";
                }
            }
            else {
                for (const bod_t &bod : bods) { // simple case of all bodies lasting from beginning to end
                    out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                           "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,"
                           "kinetic_energy,potential_energy,total_energy\r\n";
                    for (i = 0; i <= num_max_els; ++i)
                        out << i*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                            << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                            << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                            << energy[bod.id][i] << "\r\n";
                    out << ",,,,,,,,,\r\n";
                }
            }
            out << "time_elapsed,system_PE,system_KE,system_E\r\n";
            for (i = 0; i <= num_max_els; ++i)
                out << i*prev_dt*memFreq << ',' << tot_pe[i] << ',' << tot_ke[i] << ',' << tot_e[i] << "\r\n";
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
        ~system() {
            if constexpr (collisions) {
                if constexpr (memFreq)
                    delete del_bods_m;
                else
                    delete del_bods_n;
            }
            delete [] G_vals;
        }
        const bod_t &operator[](vec_size_t index) const {
            if (index >= bods.size())
                throw std::out_of_range("The requested body does not exist (index out of range).\n");
            return bods[index];
        }
        // implicitly-declared copy assignment operator is deleted, so must define my own:
        sys_t &operator=(const sys_t &other) {
            if (this == &other)
                return *this;
            this->dt = other.dt;
            this->half_dt = other.half_dt;
            this->iterations = other.iterations;
            this->G = other.G;
            this->bods = other.bods;
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                this->clear_evolution();
            return *this;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t,
                  bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        sys_t &operator=(const system<m, r, t, prg, mrg, coll, mF, fF, bF> &other) {
            if constexpr (std::same_as<sys_t, decltype(other)>) {
                if (this == &other)
                    return *this;
            }
            this->dt = other.dt;
            this->half_dt = other.half_dt;
            this->iterations = other.iterations;
            this->G = other.G;
            this->bods.assign(other.bods.begin(), other.bods.end());
            copy(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq) {
                this->clear_evolution();
                // if constexpr (mF) {
                //     this->del_bods_m->insert(other.del_bods_m.begin(), other.del_bods_m.end());
                // }
            }
            return *this;
        }
        sys_t &operator=(sys_t &&other) noexcept {
            if (this == &other)
                return *this;
            this->dt = std::move(other.dt);
            this->half_dt = other.half_dt;
            this->iterations = std::move(other.iterations);
            this->G = std::move(other.G);
            this->bods.assign(std::make_move_iterator(other.bods.begin()), std::make_move_iterator(other.bods.end()));
            move(this->G_vals, other.G_vals, 6*sizeof(uint64_t));
            if constexpr (memFreq)
                this->clear_evolution();
            return *this;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t,
                  bool prg, bool mrg, int coll, uint64_t mF, uint64_t fF, bool bF>
        friend std::ostream &operator<<(std::ostream&, const system<m, r, t, prg, mrg, coll, mF, fF, bF>&);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg1, bool prg2, bool mrg1, bool mrg2,
                  int c1, int c2, uint64_t mF1, uint64_t mF2, uint64_t fF1, uint64_t fF2, bool bF1, bool bF2>
        friend system<m, r, t, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
        operator+(const system<m, r, t, prg1, mrg1, c1, mF1, fF1, bF1>&,
                  const system<m, r, t, prg2, mrg2, c2, mF2, fF2, bF2>&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, uint64_t, uint64_t, bool>
        friend class system;
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                isNumWrapper, bool, uint64_t>
        friend class astro_scene;
    }; // all these template parameters are driving me coocoo
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T,bool prg,bool mrg,int coll,uint64_t mF,uint64_t fF,bool bF>
    std::ostream &operator<<(std::ostream &os, const system<M, R, T, prg, mrg, coll, mF, fF, bF> &sys) {
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
             uint64_t mF1, uint64_t mF2, uint64_t fF1, uint64_t fF2, bool bF1, bool bF2>
    system<M, R, T, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
    operator+(const system<M, R, T, prg1, mrg1, c1, mF1, fF1, bF1> &sys1,
              const system<M, R, T, prg2, mrg2, c2, mF2, fF2, bF2> &sys2) {
        typename std::vector<body<M, R, T, mF1>>::size_type size1 = sys1.bods.size();
        typename std::vector<body<M, R, T, mF2>>::size_type size2 = sys2.bods.size();
        if ((!size1 && !size2) || sys1.G != sys2.G)
            return {{}};
        if (!size1)
            return {sys2};
        if (!size2)
            return {sys1};
        system<M, R, T, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
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
        copy(ret_sys.G_vals, sys1.G_vals, 6*sizeof(uint64_t));
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
    typedef system<long double, long double, long double, false, false, 0, 0, 0, false> sys_n;
    typedef system<long double, long double, long double, true, false, 0, 0, 0, false> sys_p_n;
    typedef system<long double, long double, long double, false, true, 0, 0, 0, false> sys_m_n;
    typedef system<long double, long double, long double, false, false, 3, 0, 0, false> sys_c_n;
    typedef system<long double, long double, long double, false, false, 7, 0, 0, false> sys_C_n;
    typedef system<long double, long double, long double, true, true, 0, 0, 0, false> sys_pm_n;
    typedef system<long double, long double, long double, true, false, 3, 0, 0, false> sys_pc_n;
    typedef system<long double, long double, long double, true, false, 7, 0, 0, false> sys_pC_n;
    typedef system<long double, long double, long double, false, true, 3, 0, 0, false> sys_mc_n;
    typedef system<long double, long double, long double, false, true, 7, 0, 0, false> sys_mC_n;
    typedef system<long double, long double, long double, true, true, 3, 0, 0, false> sys_pmc_n;
    typedef system<long double, long double, long double, true, true, 7, 0, 0, false> sys_pmC_n;
    typedef system<long double, long double, long double, false, false, 0, 0, 1, false> sys_txt;
    typedef system<long double, long double, long double, true, false, 0, 0, 1, false> sys_p_txt;
    typedef system<long double, long double, long double, false, true, 0, 0, 1, false> sys_m_txt;
    typedef system<long double, long double, long double, false, false, 3, 0, 1, false> sys_c_txt;
    typedef system<long double, long double, long double, false, false, 7, 0, 1, false> sys_C_txt;
    typedef system<long double, long double, long double, true, true, 0, 0, 1, false> sys_pm_txt;
    typedef system<long double, long double, long double, true, false, 3, 0, 1, false> sys_pc_txt;
    typedef system<long double, long double, long double, true, false, 7, 0, 1, false> sys_pC_txt;
    typedef system<long double, long double, long double, false, true, 3, 0, 1, false> sys_mc_txt;
    typedef system<long double, long double, long double, false, true, 7, 0, 1, false> sys_mC_txt;
    typedef system<long double, long double, long double, true, true, 3, 0, 1, false> sys_pmc_txt;
    typedef system<long double, long double, long double, true, true, 7, 0, 1, false> sys_pmC_txt;
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
    typedef system<long double, long double, long double, false, false, 0, 1, 1, false> sys_MF_txt;
    typedef system<long double, long double, long double, true, false, 0, 1, 1, false> sys_p_MF_txt;
    typedef system<long double, long double, long double, false, true, 0, 1, 1, false> sys_m_MF_txt;
    typedef system<long double, long double, long double, false, false, 3, 1, 1, false> sys_c_MF_txt;
    typedef system<long double, long double, long double, false, false, 7, 1, 1, false> sys_C_MF_txt;
    typedef system<long double, long double, long double, true, true, 0, 1, 1, false> sys_pm_MF_txt;
    typedef system<long double, long double, long double, true, false, 3, 1, 1, false> sys_pc_MF_txt;
    typedef system<long double, long double, long double, true, false, 7, 1, 1, false> sys_pC_MF_txt;
    typedef system<long double, long double, long double, false, true, 3, 1, 1, false> sys_mc_MF_txt;
    typedef system<long double, long double, long double, false, true, 7, 1, 1, false> sys_mC_MF_txt;
    typedef system<long double, long double, long double, true, true, 3, 1, 1, false> sys_pmc_MF_txt;
    typedef system<long double, long double, long double, true, true, 7, 1, 1, false> sys_pmC_MF_txt;
}
#undef MEM_LOOP
#undef FUNC_TEMPL_SELECT
#undef EMPTY
#endif
