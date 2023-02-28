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
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              bool prog = false, bool mergeOverlappingBodies = false, int collisions = 0,
              ull_t memFreq = 0, ull_t fileFreq = 1, bool binaryFile = true>
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
        static constexpr int barnes_hut = 65536; // static constants to determine the force approx. method to use
        static constexpr int fast_multipole = 131072;
        static constexpr int no_coll_check = 0;
        static constexpr int overlap_coll_check = 3;
        static constexpr int pred_coll_check = 7;
        static constexpr long double G_SI = 0.00000000006743; // G in SI units (m^3 kg^-1 s^-2)
    private: // G, below, has not been made a static constant as it varies between objects, depending on units passed
        long double G = 66743; // Newtonian constant of Gravitation (10^(-15) m^3 kg^-1 s^-2)
        using bod_t = body<M, R, T, memFreq>;
        using vec_size_t = typename std::vector<bod_t>::size_type;
        using sys_t = system<M, R, T, prog, mergeOverlappingBodies, collisions, memFreq, fileFreq, binaryFile>;
        using mom_t = decltype(M{}*T{});
        std::vector<bod_t> bods; // not using set or map as need fast random access to bodies
        std::set<bod_t, std::less<>> del_bods; // set to store bodies removed from system after collision mergers
        long double dt = 1;
        long double half_dt = dt/2; // I gave it its own variable, since it is a commonly used quantity
        ull_t iterations = 1000;
        long double prev_dt{}; // used in methods called after evolve(), since could be changed by setter
        long double time_elapsed{};
        ull_t prev_iterations{}; // same here
        mom_t min_tot_com_mom{BILLION*10}; // minimum sum of magnitudes of COM momenta for two bodies to merge
        // using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        ull_t steps{}; // defined as an instance variable since it's required in numerous functions
        T total_pe{}; // these 5 variables are defined as instance variables to avoid their redefinition in many funcs
        T total_ke{};
        ull_t inner{};
        ull_t outer{};
        vec_size_t num_bods{};
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
        static inline ull_t to_ull(String &&str) {
            if (!str.isnumeric())
                return -1;
            // str.strip();
            ull_t total = 0;
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
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    if constexpr (mergeOverlappingBodies)
                        if (inner == outer) // will be true at some point if there has been a merger
                            continue;
                    inner_b = bods.data() + inner;
                    if ((outer_b->curr_pos - inner_b->curr_pos).magnitude() < outer_b->radius + inner_b->radius) {
                        if constexpr (!mergeOverlappingBodies) {
                            String str = "The bodies with id=";
                            str.append_back(outer_b->id).append_back(" and id=").append_back(inner_b->id);
                            str.append_back(" that were added to this system object overlap.\n");
                            throw overlapping_bodies_error(str.c_str());
                        }
                        *outer_b += *inner_b; // merges the two overlapping bodies
                        bods.erase(bods.begin() + inner); // thus the number of total bodies is reduced by 1
                        inner = 0; // have to recalculate possible mergers for newly created body
                        --num_bods;
                    }
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
            unsigned long long count = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += half_dt*(bod.curr_vel + std::get<1>(predicted_vals[count]));
                bod.curr_vel += half_dt*(bod.acc + std::get<2>(predicted_vals[count++]));
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
            unsigned long long count = 0;
            for (bod_t &bod : bods) {
                bod.curr_pos += dt*std::get<1>(predicted_vals[count]);
                bod.curr_vel += dt*std::get<2>(predicted_vals[count++]);
                if constexpr (mem)
                    bod.add_pos_vel_ke();
                if constexpr (file) {
                    // WRITE TO FILE
                }
            }
        }
        void create_energy_vectors() requires (memFreq != 0) {
            unsigned long long iters_p1 = (iterations + 1)/memFreq;
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
            for (; outer < num_bods; ++outer) {
                merging_bod = nullptr;
                overlapping.clear();
                bod_o = bods.data() + outer;
                for (inner = outer + 1; inner < num_bods; ++inner) {
                    bod_i = bods.data() + inner;
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
                    bod_t &&merged = *bod_o + *merging_bod; // create body that is the result of the merger
                    if constexpr (mem) {
                        merged.add_pos_vel_ke();
                    }
                    if constexpr (file) {
                        // WRITE TO FILE
                    }
                    if constexpr (memFreq) {
                        pe.emplace(merged.id, std::vector<T>{}); // add a std::vector to store its potential energies
                        energy.emplace(merged.id, std::vector<T>{}); // same for its total energy (it stores its own KE)
                    }
                    del_bods.emplace(std::move(*bod_o)); // move the merged bodies into the deleted bodies std::set
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
                if constexpr (mem) {
                    bod_o->add_pos_vel_ke();
                }
                if constexpr (file) {
                    // WRITE TO FILE
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
            static_assert(collisions <= pred_coll_check && !(collisions & (collisions + 1)),
                          "Invalid collision-checking option.");
        }
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
        constexpr system() : G{G_SI} {this->check_coll();}
        system(const std::initializer_list<bod_t> &list) : bods{list}, G{G_SI} {
            check_overlap();
            check_coll();
            if constexpr (memFreq)
                clear_bodies();
        }
        system(std::initializer_list<bod_t> &&list) : bods{std::move(list)}, G{G_SI} {
            check_overlap();
            check_coll();
            if constexpr (memFreq)
                clear_bodies();
        }
        explicit system(long double timestep, ull_t num_iterations, const char *units_format) :
                dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            check_coll();
            parse_units_format(units_format);
        }
        template <ull_t mF>
        system(const std::vector<body<M, R, T, mF>> &bodies, long double timestep = 1,
               ull_t num_iterations = 1000, const char *units_format = "M1:1,D1:1,T1:1") :
               bods{bodies.begin(), bodies.end()}, dt{timestep}, iterations{num_iterations} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            check_overlap();
            check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq && mF) // no need to clear bodies if the ones passed do not record history
                clear_bodies();
        }
        system(std::vector<bod_t> &&bodies, long double timestep = 1,
               ull_t num_iterations = 1000, const char *units_format = "M1:1,D1:1,T1:1") :
               bods{std::move(bodies)}, dt{timestep}, iterations{num_iterations} {
            check_overlap();
            check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq)
                clear_bodies();
        }
        template <ull_t mF>
        system(std::vector<body<M, R, T, mF>> &&bodies, long double timestep = 1,
               ull_t num_iterations = 1000, const char *units_format = "M1:1,D1:1,T1:1") :
               dt{timestep}, iterations{num_iterations} {
            std::for_each(bodies.begin(), bodies.end(),
                          [this](body<M, R, T, mF> &b){bods.emplace_back(std::move(b));});
            check_overlap();
            check_coll();
            parse_units_format(units_format);
            if constexpr (memFreq && mF)
                clear_bodies();
        }
        /* copy constructors are made to only copy the variables seen below: it does not make sense for a new system
         * object (even when copy-constructed) to be "evolved" if it has not gone through the evolution itself */
        template <bool prg, bool mrg, int coll, ull_t mF, ull_t fF, bool bF>
        // so a system can be constructed from another without same checks
        system(const system<M, R, T, prg, mrg, coll, mF, fF, bF> &other) :
               bods{other.bods.begin(), other.bods.end()}, dt{other.dt}, iterations{other.iterations}, G{other.G} {
            check_overlap();
            check_coll();
            if constexpr (memFreq && mF)
                clear_bodies();
        }
        template <bool prg, bool mrg, int coll, ull_t fF, bool bF>
        system(system<M, R, T, prg, mrg, coll, memFreq, fF, bF> &&other) :
               bods{std::move(other.bods)}, dt{other.dt}, iterations{other.iterations}, G{other.G} {
            check_overlap();
            check_coll();
            if constexpr (memFreq)
                clear_bodies();
        }
        template <bool prg, bool mrg, int coll, ull_t mF, ull_t fF, bool bF>
        system(system<M, R, T, prg, mrg, coll, mF, fF, bF> &&other) :
               dt{other.dt}, iterations{other.iterations}, G{other.G} {
            std::for_each(other.bods.begin(), other.bods.end(),
                          [this](body<M, R, T, mF> &b){bods.emplace_back(std::move(b));});
            check_overlap();
            check_coll();
            if constexpr (memFreq && mF)
                clear_bodies();
        }
        vec_size_t num_bodies() {
            return bods.size();
        }
        bool set_iterations(const ull_t &number) noexcept {
            if (!number) // cannot have zero iterations
                return false;
            iterations = number;
            return true;
        }
        bool set_iterations(const ull_t &&number) noexcept {
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
        const bod_t &get_body(unsigned long long id) const {
            for (const auto &b : *this) // this is where it would be nice to be using a set!!
                if (b.id == id)
                    return b;
            throw std::invalid_argument("The id passed does not correspond to any body present in this system "
                                        "object.\n");
        }
        bool remove_body(unsigned long long id) {
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
        void clear_evolution() requires (memFreq != 0) {
            // for (bod_t &bod : bods)
            //     bod.clear();
            std::for_each(bods.begin(), bods.end(), [](bod_t &bod){bod.clear();});
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
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
            ull_t i;
            const ull_t num_max_els = (prev_iterations + 1)/memFreq; // CHECK THIS IS CORRECT
            if constexpr (collisions) {

            }
            for (const bod_t &bod : del_bods) {
                out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                    << "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy,"
                       "potential_energy,total_energy\r\n";
                for (i = 0; i <= num_max_els; ++i)
                    out << i*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                        << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                        << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                        << energy[bod.id][i] << "\r\n";
                out << ",,,,,,,,,\r\n";
                // ++count;
            }
            for (const bod_t &bod : bods) {
                out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                    << "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy,"
                       "potential_energy,total_energy\r\n";
                for (i = 0; i <= num_max_els; ++i)
                    out << i*prev_dt*memFreq << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                        << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                        << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe[bod.id][i] << ','
                        << energy[bod.id][i] << "\r\n";
                out << ",,,,,,,,,\r\n";
                // ++count;
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
        const bod_t &operator[](vec_size_t index) const {
            if (index >= bods.size())
                throw std::out_of_range("The requested body does not exist (index out of range).\n");
            return bods[index];
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t,
                  bool prg, bool mrg, int coll, ull_t mF, ull_t fF, bool bF>
        sys_t &operator=(const system<m, r, t, prg, mrg, coll, mF, fF, bF> &other) {
            if (this == &other)
                return *this;
            this->dt = other.dt;
            this->iterations = other.iterations;
            this->G = other.G;
            this->bods = other.bods;
            if constexpr (memFreq)
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
            if constexpr (memFreq)
                this->clear_evolution();
            return *this;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t,
                  bool prg, bool mrg, int coll, ull_t mF, ull_t fF, bool bF>
        friend std::ostream &operator<<(std::ostream&, const system<m, r, t, prg, mrg, coll, mF, fF, bF>&);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg1, bool prg2, bool mrg1, bool mrg2,
                  int c1, int c2, ull_t mF1, ull_t mF2, ull_t fF1, ull_t fF2, bool bF1, bool bF2>
        friend system<m, r, t, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
        operator+(const system<m, r, t, prg1, mrg1, c1, mF1, fF1, bF1>&,
                  const system<m, r, t, prg2, mrg2, c2, mF2, fF2, bF2>&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, ull_t, ull_t, bool>
        friend class system;
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                isNumWrapper, bool, ull_t>
        friend class astro_scene;
    }; // all these template parameters are driving me coocoo
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, bool prg, bool mrg, int coll, ull_t mF, ull_t fF, bool bF>
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
             ull_t mF1, ull_t mF2, ull_t fF1, ull_t fF2, bool bF1, bool bF2>
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
