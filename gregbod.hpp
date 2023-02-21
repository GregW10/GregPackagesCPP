//
// Created by mario on 08/10/2022.
//

#ifndef GREGBOD_H
#define GREGBOD_H

#include "gregvec.h"
#include <set>
#include <fstream>
#include <map>
#include <regex>

namespace gtd {
    class nbody_error : public std::logic_error {
    public:
        nbody_error() : std::logic_error{"Invalid parameters for an N-body simulation.\n"} {}
        explicit nbody_error(const char *message) : std::logic_error{message} {}
    };
    class overlapping_bodies_error : public nbody_error {
        const char *message;
        String str;
    public:
        overlapping_bodies_error() {
            message = "Two gtd::system objects cannot be added that will produce an overlap of 2 or more bodies "
                      "within.\n";
        }
        overlapping_bodies_error(unsigned long long id1, unsigned long long id2) {
            str.append_back("The addition of these two gtd::system objects would produce an overlap between the body "
                            "with id=").append_back(id1).append_back(" and the body with id=").append_back(id2).
                            append_back(".\nThe addition of two gtd::system objects cannot result in overlap between "
                                        "any bodies.\n");
            message = str.c_str();
        }
        explicit overlapping_bodies_error(const char *msg) : message(msg) {}
        const char *what() const noexcept override {
            return message;
        }
    };
    class negative_mass_error : public nbody_error {
        const char *message;
    public:
        negative_mass_error() : message("A gtd::body cannot have a negative mass.\n") {}
        explicit negative_mass_error(const char *msg) : message(msg) {}
        const char *what() const noexcept override {
            return message;
        }
    };
    class negative_radius_error : public nbody_error {
        const char *message;
    public:
        negative_radius_error() : message("A gtd::body cannot have a negative radius.\n") {}
        explicit negative_radius_error(const char *msg) : message(msg) {}
        const char *what() const noexcept override {
            return message;
        }
    };
    class two_body_error : public nbody_error {
    public:
        two_body_error() : nbody_error{"A two-body integration can only be performed if there are two bodies present in"
                                       " the system.\n"} {}
        explicit two_body_error(const char *message) : nbody_error{message} {}
    };
    class empty_system_error : public nbody_error {
    public:
        empty_system_error() : nbody_error{"This operation is not permitted on empty system objects (no bodies).\n"} {}
        explicit empty_system_error(const char *message) : nbody_error{message} {}
    };
    class no_evolution_error : public nbody_error {
    public:
        no_evolution_error() : nbody_error{"This operation is not permitted on system objects that have not been "
                                           "evolved.\n"} {}
        explicit no_evolution_error(const char *message) : nbody_error{message} {}
    };
    class invalid_id_error : public nbody_error {
        unsigned long long invalid_id{};
        String msg;
    public:
        invalid_id_error() : nbody_error{"Body with specified ID was not found.\n"} {}
        explicit invalid_id_error(const char *message) : nbody_error{message} {}
        explicit invalid_id_error(unsigned long long id) : invalid_id{id} {
            msg.append_back("No body found with ID = ").append_back(invalid_id).append_back(".\n");
        }
        const char *what() const noexcept override {
            if (msg == "")
                return nbody_error::what();
            return msg.c_str();
        }
    };
    template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, ull_t, ull_t, bool>
    class system;
    /* body_counter was created because each instantiation of the body subclass with different template parameters is
     * actually a different class, so the count of bodies for each different class template instantiation would be
     * different */
    class body_counter {
        static inline ull_t count = 0;
        static inline std::set<ull_t> ids;
        /* the below mutex is present to allow body_counter objects to be created concurrently in different threads,
         * without two IDs ever being the same */
        static inline std::mutex id_mutex;
        // void set_id() {
        //     unsigned long long prev = 0;
        //     for (const auto &i : ids) {
        //         if (i - prev > 1) {
        //             id = prev + 1;
        //             return;
        //         }
        //     }
        //     id = *(ids.end()--) + 1;
        // }
    protected:
        ull_t id; // is immutable after the object has been constructed - unless the object is moved
    public:
        body_counter() {
            // id = count;
            // if (ids.contains(id))
            //     set_id();
            this->id = 0;
            std::lock_guard<std::mutex> lock{id_mutex};
            while (ids.contains(id)) ++this->id;
            ids.insert(this->id);
            ++count;
        }
        body_counter(body_counter &&other) noexcept {
            this->id = other.id;
            other.id = -1; // would never, ever be reached (normally)
            std::lock_guard<std::mutex> lock{id_mutex};
            while (ids.contains(other.id))
                --other.id;
            ids.insert(other.id); // important in case other is not actually destroyed after this call
            ++count; // this will prob. cause count to actually be the same as before the call to the constructor, since
        } // ... "other" is likely destroyed at the end of the call, reducing count by 1, so count must be incremented.
        body_counter(const body_counter &other) = delete; // again, all IDs must be unique
        ull_t get_id() const noexcept {
            return id;
        }
        static ull_t get_tot_bc() noexcept {
            return count;
        }
        static std::vector<ull_t> get_all_ids() {
            return {ids.begin(), ids.end()};
        }
        static void print_all_ids() {
            ull_t b_count = 0;
            for (const ull_t &i : ids)
                printf("ID: %llu, count: %llu/%llu\n", i, ++b_count, count);
        }
        ~body_counter() {
            std::lock_guard<std::mutex> lock{id_mutex};
            --count;
            ids.erase(this->id);
        }
        body_counter &operator=(const body_counter&) = delete; // all IDs must be unique
        body_counter &operator=(body_counter &&other) noexcept {
            this->id = other.id;
            other.id = -1;
            std::lock_guard<std::mutex> lock{id_mutex};
            while (ids.contains(other.id))
                --other.id;
            ids.insert(other.id);
            return *this;
        }
        explicit operator unsigned long long() const noexcept {
            return this->id; // thus, if desired, a body_counter can be cast to an unsigned long long to get its ID
        }
        /* the first three overloads (for <) are essential as they allow a transparent comparator to be used within any
         * std::set used to store bodies (so that lookup can be performed based on ID and not on a body itself, thus
         * allowing a std::set to be used instead of a std::map) */
        friend bool operator<(const body_counter&, const body_counter&);
        friend bool operator<(const unsigned long long&, const body_counter&);
        friend bool operator<(const body_counter&, const unsigned long long&);
        friend bool operator>(const body_counter&, const body_counter&);
        friend bool operator>(const unsigned long long&, const body_counter&);
        friend bool operator>(const body_counter&, const unsigned long long&);
        friend bool operator<=(const body_counter&, const body_counter&);
        friend bool operator<=(const unsigned long long&, const body_counter&);
        friend bool operator<=(const body_counter&, const unsigned long long&);
        friend bool operator>=(const body_counter&, const body_counter&);
        friend bool operator>=(const unsigned long long&, const body_counter&);
        friend bool operator>=(const body_counter&, const unsigned long long&);
        friend bool operator==(const body_counter&, const body_counter&);
        friend bool operator==(const unsigned long long&, const body_counter&);
        friend bool operator==(const body_counter&, const unsigned long long&);
        friend bool operator!=(const body_counter&, const body_counter&);
        friend bool operator!=(const unsigned long long&, const body_counter&);
        friend bool operator!=(const body_counter&, const unsigned long long&);
        friend std::ostream &operator<<(std::ostream&, const body_counter&);
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, ull_t, ull_t, bool>
        friend class system;
    };
    bool operator<(const body_counter &bc1, const body_counter &bc2) {
        return bc1.id < bc2.id;
    }
    bool operator<(const unsigned long long &ID, const body_counter &bc) {
        return ID < bc.id;
    }
    bool operator<(const body_counter &bc, const unsigned long long &ID) {
        return bc.id < ID;
    }
    bool operator>(const body_counter &bc1, const body_counter &bc2) {
        return bc1.id > bc2.id;
    }
    bool operator>(const unsigned long long &ID, const body_counter &bc) {
        return ID > bc.id;
    }
    bool operator>(const body_counter &bc, const unsigned long long &ID) {
        return bc.id > ID;
    }
    bool operator<=(const body_counter &bc1, const body_counter &bc2) {
        return bc1.id <= bc2.id;
    }
    bool operator<=(const unsigned long long &ID, const body_counter &bc) {
        return ID <= bc.id;
    }
    bool operator<=(const body_counter &bc, const unsigned long long &ID) {
        return bc.id <= ID;
    }
    bool operator>=(const body_counter &bc1, const body_counter &bc2) {
        return bc1.id >= bc2.id;
    }
    bool operator>=(const unsigned long long &ID, const body_counter &bc) {
        return ID >= bc.id;
    }
    bool operator>=(const body_counter &bc, const unsigned long long &ID) {
        return bc.id >= ID;
    }
    bool operator==(const body_counter &bc1, const body_counter &bc2) {
        return bc1.id == bc2.id;
    }
    bool operator==(const unsigned long long &ID, const body_counter &bc) {
        return ID == bc.id;
    }
    bool operator==(const body_counter &bc, const unsigned long long &ID) {
        return bc.id == ID;
    }
    bool operator!=(const body_counter &bc1, const body_counter &bc2) {
        return bc1.id != bc2.id;
    }
    bool operator!=(const unsigned long long &ID, const body_counter &bc) {
        return ID != bc.id;
    }
    bool operator!=(const body_counter &bc, const unsigned long long &ID) {
        return bc.id != ID;
    }
    std::ostream &operator<<(std::ostream &os, const body_counter &bc) {
        return os << "[gtd::body_counter@" << &bc << ":id=" << bc.id << ']';
    }
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double,
              unsigned long long recFreq = 1>
    class body : public body_counter {
        /* A body object is one which represents a real, spherical object in 3D space. It has a given mass & radius, &
         * position & velocity at any given time. It is able to record its entire history of where it has been. */
        /* The template parameters M, R & T are the data types that are used to represent the body's mass, radius &
         * position & velocity, respectively. */
        /* The non-type template parameter recFreq represents the frequency with which a body should record its
         * position, velocity & kinetic energy. In other words, it dictates how many shifts/updates in position &
         * velocity take place per single recording. A value of zero indicates no recording. */
        /* A body is not thread-safe and should never be modified by multiple threads. */
    public:
        static inline bool allow_self_addition = false;
    protected:
        R radius{1};
        vector3D<T> curr_pos{}; // current position
    private:
        M mass_{1};
        vector3D<T> curr_vel{}; // current velocity - have not included acc. as this should always be det. externally
        vector3D<T> acc{}; // current acceleration of the body - only used within system class
        T pe{}; // current potential energy of the body - only used in system class since requires other bodies to calc.
        using K = decltype(0.5*mass_*curr_vel.magnitude()*curr_vel.magnitude());
        K curr_ke{}; // current kinetic energy
        long double rest_c{1}; // coeff. of restitution - between 0-1 (0 = perfectly inelastic, 1 = perfectly elastic)
        /* Only quantities that represent the state of a body on its own are stored (below) - position, velocity, and
         * kinetic energy - whilst quantities that depend on the states of other bodies are not stored - acceleration
         * and potential energy. */
        std::vector<vector3D<T>> positions; // these std::vectors will hold all the positions and velocities of the
        std::vector<vector3D<T>> velocities; // ... body as it moves
        std::vector<K> energies;
        unsigned long long rec_counter = 0;
        mutable std::vector<std::tuple<const vector3D<T>&, const vector3D<T>&, const K&>> cref;
        typedef typename std::vector<vector3D<T>>::size_type vec_s_t;
        static inline String def_path = get_home_path<String>() + FILE_SEP + "Body_Trajectory_";
        void add_pos_vel() {
            positions.push_back(curr_pos);
            velocities.push_back(curr_vel);
        }
        void add_pos_vel_ke() {
            add_pos_vel();
            set_ke();
            energies.push_back(curr_ke);
        }
        void emplace_zeros() {
            this->positions.emplace_back();
            this->velocities.emplace_back();
            this->energies.emplace_back();
        }
        void set_ke() {
            auto &&mag = curr_vel.magnitude(); // best to push var. onto stack rather than call func. twice
            curr_ke = 0.5*mass_*mag*mag;
        }
        void check_mass() { // I prefer to throw an exception, rather than simply taking no action, to make it clear
            if (mass_ < M{0}) // ... where a negative quantity has attempted to be set
                throw negative_mass_error();
        }
        void check_radius() const {
            if (radius < R{0})
                throw negative_radius_error();
        }
        body(M &&body_mass, R &&body_radius, vector3D<T> &&pos, vector3D<T> &&vel, vector3D<T> &&acceleration,
             const long double &&restitution) :
        mass_{std::move(body_mass)}, radius{std::move(body_radius)}, curr_pos{std::move(pos)}, curr_vel{std::move(vel)},
        acc{acceleration}, rest_c{restitution}
        {add_pos_vel_ke(); check_mass(); check_radius();} // private ctor as only needed for 1 func.
    public:
        /* unfortunately, the constructors cannot be marked constexpr, because it is impossible for any body_counter
         * constructor to be constexpr (since ID must be dynamically determined) */
        body() {
            if constexpr (recFreq)
                this->emplace_zeros();
        }
        body(M &&body_mass, R &&body_radius) : mass_{std::move(body_mass)}, radius{std::move(body_radius)} {
            check_mass();
            check_radius();
            if constexpr (recFreq)
                this->emplace_zeros();
        }
        body(const M &body_mass, const R &body_radius) : mass_{body_mass}, radius{body_radius} {
            check_mass();
            check_radius();
            if constexpr (recFreq)
                this->emplace_zeros();
        }
        body(M &&body_mass, R &&body_radius, vector3D<T> &&pos, vector3D<T> &&vel) :
        mass_{std::move(body_mass)}, radius{std::move(body_radius)}, curr_pos{std::move(pos)}, curr_vel{std::move(vel)}{
            check_mass();
            check_radius();
            if constexpr (recFreq)
                this->add_pos_vel_ke();
        }
        body(const M &body_mass, const R &body_radius, const vector3D<T> &pos, const vector3D<T> &vel) :
        mass_{body_mass}, radius{body_radius}, curr_pos{pos}, curr_vel{vel} {
            check_mass();
            check_radius();
            if constexpr (recFreq)
                this->add_pos_vel_ke();
        }
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t, ull_t rF>
        body(const body<m, r, t, rF> &other) :
        mass_{other.mass_}, radius{other.radius}, curr_pos{other.curr_pos}, curr_vel{other.curr_vel},
        curr_ke{other.curr_ke}, acc{other.acc}, rest_c{other.rest_c}, pe{other.pe} {
            if constexpr (recFreq) {
                this->add_pos_vel();
                energies.push_back(this->curr_ke);
            }
        }
        template <ull_t rF>
        body(const body<M, R, T, rF> &other) :
        mass_{other.mass_}, radius{other.radius}, curr_pos{other.curr_pos},
        curr_vel{other.curr_vel}, curr_ke{other.curr_ke}, acc{other.acc}, rest_c{other.rest_c}, pe{other.pe} {
            if constexpr (recFreq) {
                this->add_pos_vel();
                energies.push_back(this->curr_ke);
            }
        }
        body(body<M, R, T, recFreq> &&other) noexcept :
             body_counter{std::move(other)}, mass_{std::move(other.mass_)},
             radius{std::move(other.radius)}, curr_pos{std::move(other.curr_pos)}, curr_vel{std::move(other.curr_vel)},
             curr_ke{std::move(other.curr_ke)}, acc{std::move(other.acc)}, positions{std::move(other.positions)},
             velocities{std::move(other.velocities)}, energies{std::move(other.energies)}, cref{std::move(other.cref)},
        rest_c{other.rest_c}, pe{other.pe}, rec_counter{other.rec_counter} {
            if constexpr (recFreq)
                if (rec_counter == recFreq) {
                    this->add_pos_vel_ke();
                    rec_counter = 0;
                }
        }
        template <ull_t rF>
        body(body<M, R, T, rF> &&other) noexcept :
             body_counter{std::move(other)}, mass_{std::move(other.mass_)},
             radius{std::move(other.radius)}, curr_pos{std::move(other.curr_pos)}, curr_vel{std::move(other.curr_vel)},
             curr_ke{std::move(other.curr_ke)}, acc{std::move(other.acc)}, rest_c{other.rest_c}, pe{other.pe} {
            if constexpr (recFreq)
                this->add_pos_vel_ke();
        }
        const M &mass() const noexcept {
            return mass_;
        }
        const R &rad() const noexcept {
            return radius;
        }
        const vector3D<T> &pos() const noexcept {
            return curr_pos;
        }
        const vector3D<T> &vel() const noexcept {
            return curr_vel;
        }
        vector3D<T> &acceleration() noexcept {
            return acc;
        }
        const vector3D<T> &acceleration() const noexcept {
            return acc;
        }
        const K &ke() const noexcept {
            return curr_ke;
        }
        const T &potential_energy() const noexcept {
            return pe;
        }
        long double restitution() const noexcept {
            return rest_c;
        }
        const vector3D<T> &prev_pos_at(vec_s_t index) const requires (recFreq != 0) {
            // static_assert(recFreq, "prev_pos_at cannot be called on body objects with recFreq == 0 (no history)\n");
            if (index >= positions.size())
                throw std::out_of_range("Requested position does not exist (index out of range).\n");
            return positions[index];
        }
        const vector3D<T> &prev_vel_at(vec_s_t index) const requires (recFreq != 0) {
            // static_assert(recFreq, "prev_vel_at cannot be called on body objects with recFreq == 0 (no history)\n");
            if (index >= velocities.size())
                throw std::out_of_range("Requested velocity does not exist (index out of range).\n");
            return velocities[index];
        }
        const K &prev_ke_at(typename std::vector<K>::size_type index) const requires (recFreq != 0) {
            // static_assert(recFreq, "prev_ke_at cannot be called on body objects with recFreq == 0 (no history)\n");
            if (index >= energies.size())
                throw std::out_of_range("Requested kinetic energy does not exist (index out of range).\n");
            return energies[index];
        }
        bool set_mass(const M &new_mass) {
            if (new_mass < 0)
                return false;
            mass_ = new_mass;
            return true;
        }
        bool set_mass(M &&new_mass) noexcept {
            if (new_mass < 0)
                return false;
            mass_ = std::move(new_mass);
            return true;
        }
        bool set_radius(const R &new_radius) {
            if (new_radius < 0)
                return false;
            radius = new_radius;
            return true;
        }
        bool set_radius(R &&new_radius) noexcept {
            if (new_radius < 0)
                return false;
            radius = std::move(new_radius);
            return true;
        }
        void set_acc(const vector3D<T> &new_acc) {
            acc = new_acc;
        }
        void set_acc(vector3D<T> &&new_acc) noexcept {
            acc = std::move(new_acc);
        }
        bool set_pe(const T &new_pe) {
            if (new_pe > 0) // gravitational PE can only be zero (at infinity) or negative
                return false;
            pe = new_pe;
            return true;
        }
        bool set_pe(T &&new_pe) noexcept {
            if (new_pe > 0)
                return false;
            pe = std::move(new_pe);
            return true;
        }
        bool set_restitution(const long double &new_restitution_coefficient) noexcept {
            if (new_restitution_coefficient < 0 || new_restitution_coefficient > 1)
                return false;
            rest_c = new_restitution_coefficient;
            return true;
        }
        bool set_restitution(const long double &&new_restitution_coefficient) noexcept {
            if (new_restitution_coefficient < 0 || new_restitution_coefficient > 1)
                return false;
            rest_c = new_restitution_coefficient;
            return true;
        }
        void update(const vector3D<T> &new_position, const vector3D<T> &new_velocity) noexcept {
            curr_pos = new_position; // no checks to be performed, these new vecs can be anything
            curr_vel = new_velocity;
            if constexpr (recFreq) {
                if (++rec_counter != recFreq)
                    return;
                rec_counter = 0;
                this->add_pos_vel_ke();
            }
        }
        void shift(const vector3D<T> &position_shift, const vector3D<T> &velocity_shift) noexcept {
            curr_pos += position_shift;
            curr_vel += velocity_shift;
            if constexpr (recFreq) {
                if (++rec_counter != recFreq)
                    return;
                rec_counter = 0;
                this->add_pos_vel_ke();
            }
        }
        void update(vector3D<T> &&new_position, vector3D<T> &&new_velocity) noexcept {
            curr_pos = std::move(new_position);
            curr_vel = std::move(new_velocity);
            if constexpr (recFreq) {
                if (++rec_counter != recFreq)
                    return;
                rec_counter = 0;
                this->add_pos_vel_ke();
            }
        }
        void shift(const vector3D<T> &&position_shift, const vector3D<T> &&velocity_shift) noexcept {
            shift(position_shift, velocity_shift);
        }
        void apply_pos_transform(const matrix<T> &transform) { // will throw if not 3x3 matrix
            curr_pos.apply(transform);
        }
        void apply_pos_transform(const matrix<T> &&transform) {
            curr_pos.apply(transform);
        }
        void apply_vel_transform(const matrix<T> &transform) {
            curr_vel.apply(transform);
        }
        void apply_vel_transform(const matrix<T> &&transform) {
            curr_vel.apply(transform);
        }
        void apply_acc_transform(const matrix<T> &transform) {
            acc.apply(transform);
        }
        void apply_acc_transform(const matrix<T> &&transform) {
            acc.apply(transform);
        }
        auto momentum() const {
            return mass_*curr_vel;
        }
        void reset_acc() {
            acc = T{0};
        } // resets the body to initial setup and clears trajectory if specified:
        void reset(bool clear_trajectory = true) requires (recFreq != 0) {
            // static_assert(recFreq, "reset cannot be called on body objects with recFrec == 0 (no history)\n");
            curr_pos = positions.front(); // if recFreq is not zero, there is guaranteed to be at least 1 element
            curr_vel = velocities.front();
            curr_ke = energies.front();
            if (clear_trajectory) {
                clear();
                return;
            }
            if (++rec_counter != recFreq)
                return;
            rec_counter = 0;
            this->add_pos_vel();
            energies.push_back(curr_ke);
        }
        void clear() requires (recFreq != 0) {
            // static_assert(recFreq, "clear cannot be called on body objects with recFrec == 0 (no history)\n");
            positions.clear();
            velocities.clear();
            energies.clear();
            cref.clear();
            this->add_pos_vel(); // clear is the only method which will overrule recFreq (values are always added)
            energies.push_back(curr_ke); // no need to calculate ke again
            rec_counter = 0;
        }
        std::vector<vector3D<T>> get_positions_cpy() const requires (recFreq != 0) {
            return {positions};
        }
        std::vector<vector3D<T>> get_velocities_cpy() const requires (recFreq != 0) {
            return {velocities};
        }
        std::vector<K> get_kinetic_energies_cpy() const requires (recFreq != 0) {
            return {energies};
        }
        const std::vector<vector3D<T>> &get_positions() const requires (recFreq != 0) {
            return positions;
        }
        const std::vector<vector3D<T>> &get_velocities() const requires (recFreq != 0) {
            return velocities;
        }
        const std::vector<K> &get_kinetic_energies() const requires (recFreq != 0) {
            return energies;
        }
        /* since one should not modify any past pos, vel or ke at all, and no current pos, vel or ke other than
         * through the update() or shift() methods, all the iteration methods below return const iterators */
        auto begin() const requires (recFreq != 0) {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).cbegin();
            return cref.cbegin();
        }
        auto end() const requires (recFreq != 0) {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).cend();
            return cref.cend();
        }
        auto cbegin() const requires (recFreq != 0) {
            return begin();
        }
        auto cend() const requires (recFreq != 0) {
            return end();
        }
        auto rbegin() const requires (recFreq != 0) {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).crbegin();
            return cref.crbegin();
        }
        auto rend() const requires (recFreq != 0) {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).crend();
            return cref.crend();
        }
        auto crbegin() const requires (recFreq != 0) {
            return rbegin();
        }
        auto crend() const requires (recFreq != 0) {
            return rend();
        }
        bool trajectory_to_txt(std::ofstream &out, bool full_csv_style = true) const requires (recFreq != 0) {
            if (!out.good())
                return false;
            ull_t count = 0;
            if (full_csv_style) {
                out << "body_id,mass,radius\r\n" << body_counter::id << ',' << mass_ << ',' << radius << "\r\n"
                    << "iteration,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy\r\n";
                for (const auto &[pos, vel, ke] : *this) {
                    out << count << ',' << pos.x << ',' << pos.y << ',' << pos.z << ',' << vel.x << ',' << vel.y
                        << ',' << vel.z << ',' << ke << "\r\n";
                    count += recFreq;
                }
                out << "\r\n";
                out.close();
                return true;
            }
            long long before = out.tellp();
            out << "Body ID: " << body_counter::id << ", Mass = " << mass_ << ", Radius = " << radius << "\n";
            long long after = out.tellp() - before;
            --after;
            for (long i = 0; i < after; ++i)
                out.put('-');
            out.put('\n');
            out << "iteration,position,velocity,kinetic_energy\n";
            for (const auto &[pos, vel, ke] : *this) {
                out << count << ',' << pos << ',' << vel << ',' << ke << '\n';
                count += recFreq;
            }
            out << '\n';
            return true;
        }
        bool trajectory_to_txt(const String &path = def_path, bool truncate = false, bool full_csv_style = false)
        const requires (recFreq != 0) {
            if (&path == &def_path)
                def_path.append_back(get_date_and_time()).append_back(".csv");
            std::ofstream out(path.c_str(), truncate ? std::ios_base::trunc : std::ios_base::app);
            bool ret = trajectory_to_txt(out, full_csv_style);
            out.close();
            if (&path == &def_path)
                def_path.erase_chars(def_path.get_length() - 29);
            return ret;
        }
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t, ull_t rF>
        body<M, R, T, recFreq> &operator=(const body<m, r, t, rF> &other) {
            if (&other == this)
                return *this;
            this->mass_ = other.mass_;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
            this->acc = other.acc;
            this->pe = other.pe;
            this->rest_c = other.rest_c;
            this->rec_counter = 0; // initial values always added (done in clear)
            clear();
            return *this;
        }
        template <ull_t rF>
        body<M, R, T, recFreq> &operator=(const body<M, R, T, rF> &other) {
            if (&other == this)
                return *this;
            this->mass_ = other.mass_;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
            this->acc = other.acc;
            this->pe = other.pe;
            this->rest_c = other.rest_c;
            this->rec_counter = 0;
            clear();
            return *this;
        }
        template <unsigned long long rF>
        body<M, R, T, recFreq> &operator=(body<M, R, T, rF> &&other) noexcept {
            if (&other == this)
                return *this;
            this->mass_ = std::move(other.mass_);
            this->radius = std::move(other.radius);
            this->curr_pos = std::move(other.curr_pos);
            this->curr_vel = std::move(other.curr_vel);
            this ->curr_ke = std::move(other.curr_ke);
            this->acc = std::move(other.acc);
            this->pe = std::move(other.pe);
            this->rest_c = other.rest_c; // is long double, so not moved
            if constexpr (recFreq) {
                this->positions = std::move(other.positions);
                this->velocities = std::move(other.velocities);
                this->energies = std::move(other.energies);
                this->cref = std::move(other.cref);
                this->rec_counter = other.rec_counter;
            }
            return *this;
        }
        vector3D<T> operator[](vec_s_t index) const requires (recFreq != 0) { // returns a copy
            if (index >= positions.size())
                throw std::out_of_range("The specified position index is out of range.\n");
            return positions[index];
        }
        template <ull_t rF>
        body<M, R, T, recFreq> &operator+=(const body<M, R, T, rF>&other)noexcept requires isConvertible<R,long double>{
            if (&other == this) {
                if (!allow_self_addition)
                    return *this;
                this->mass_ *= M{2};
                this->radius = cbrtl(2*this->radius*this->radius*this->radius);
                return *this; // clearly, for self-addition, position, velocity & acceleration remain unchanged
            }
            this->curr_pos = com(*this, other);
            this->curr_vel = vel_com(*this, other);
            this->acc = acc_com(*this, other);
            this->rest_c = MEAN_AVG(this->rest_c, other.rest_c); // average coefficient of restitution
            this->mass_ += other.mass_;
            this->radius = cbrtl(this->radius*this->radius*this->radius + other.radius*other.radius*other.radius);
            if constexpr (recFreq) {
                if (++rec_counter == recFreq) {
                    this->add_pos_vel_ke();
                    rec_counter = 0;
                }
            }
            return *this;
        }
        virtual ~body() = default;
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, ull_t rF>
        friend std::ostream &operator<<(std::ostream&, const body<m, r, t, rF>&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto com(const body<m1, r1, t1, rF1>&, const body<m2, r2, t2, rF2>&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto vel_com(const body<m1, r1, t1, rF1>&, const body<m2, r2, t2, rF2>&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto acc_com(const body<m1, r1, t1, rF1>&, const body<m2, r2, t2, rF2>&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto operator+(const body<m1, r1, t1, rF1>&, const body<m2, r2, t2, rF2>&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto operator+(const body<m1, r1, t1, rF1>&, const body<m2, r2, t2, rF2>&&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto operator+(const body<m1, r1, t1, rF1>&&, const body<m2, r2, t2, rF2>&);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
                  isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
        friend inline auto operator+(const body<m1, r1, t1, rF1>&&, const body<m2, r2, t2, rF2>&&);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, bool prg1, bool prg2, bool mrg1, bool mrg2,
                  int c1, int c2, ull_t mF1, ull_t mF2, ull_t fF1, ull_t fF2, bool bF1, bool bF2>
        friend system<m, r, t, prg1 & prg2, mrg1 & mrg2, c1 & c2, MEAN_AVG(mF1, mF2), MEAN_AVG(fF1, fF2), bF1 & bF2>
        operator+(const system<m, r, t, prg1, mrg1, c1, mF1, fF1, bF1> &sys1,
                  const system<m, r, t, prg2, mrg2, c2, mF2, fF2, bF2> &sys2);
        template <isNumWrapper, isNumWrapper, isNumWrapper>
        friend class ray;
        template <isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper, isNumWrapper,
                  isNumWrapper, bool>
        friend class astro_scene;
        template <isNumWrapper, isNumWrapper, isNumWrapper, ull_t>
        friend class body;
        template <isNumWrapper, isNumWrapper, isNumWrapper, bool, bool, int, ull_t, ull_t, bool>
        friend class system;
    };
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t, ull_t rF>
    std::ostream &operator<<(std::ostream &os, const body<m, r, t, rF> &bod) {
        return os << "[gtd::body@" << &bod << ":id=" << bod.id << ",m=" << +bod.mass_ << ",r=" << +bod.radius
                  << ",current_pos=(" << bod.curr_pos << "),current_vel=(" << bod.curr_vel << "),current_ke="
                  << bod.curr_ke << "]";
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto com(const body<m1, r1, t1, rF1> &b1, const body<m2, r2, t2, rF2> &b2) { // centre of mass_ position
        return (b1.mass_*b1.curr_pos + b2.mass_*b2.curr_pos)/(b1.mass_ + b2.mass_);
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto vel_com(const body<m1, r1, t1, rF1> &b1, const body<m2, r2, t2, rF2> &b2) {
        /* centre-of-mass velocity, taking conservation of momentum into account */
        return (b1.momentum() + b2.momentum())/(b1.mass_ + b2.mass_);
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto acc_com(const body<m1, r1, t1, rF1> &b1, const body<m2, r2, t2, rF2> &b2) {
        return (b1.mass_*b1.acc + b2.mass_*b2.acc)/(b1.mass_ + b2.mass_);
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto operator+(const body<m1, r1, t1, rF1> &b1, const body<m2, r2, t2, rF2> &b2) {
        /* performs a "merging" of two bodies, with the new body having the sum of both masses, being at the
         * centre-of-mass of both bodies, having a volume equal to the sum of both volumes (using the radii), and a new
         * velocity & acceleration such that momentum is conserved (i.e., COM velocity and acceleration) */
        /* self-addition is not considered here as the returned body is new (body has not been added onto itself) */
        return body<decltype(b1.mass_ + b2.mass_), long double, decltype(m1{}*t2{}/m2{}),
                    MEAN_AVG(rF1, rF2)>(b1.mass_ + b2.mass_,
        cbrtl(b1.radius*b1.radius*b1.radius + b2.radius*b2.radius*b2.radius), com(b1, b2), vel_com(b1, b2),
        acc_com(b1, b2), MEAN_AVG(b1.rest_c, b2.rest_c));
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto operator+(const body<m1, r1, t1, rF1> &b1, const body<m2, r2, t2, rF2> &&b2) {
        return b1 + b2;
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto operator+(const body<m1, r1, t1, rF1> &&b1, const body<m2, r2, t2, rF2> &b2) {
        return b1 + b2;
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, ull_t rF1,
              isNumWrapper m2, isNumWrapper r2, isNumWrapper t2, ull_t rF2>
    inline auto operator+(const body<m1, r1, t1, rF1> &&b1, const body<m2, r2, t2, rF2> &&b2) {
        return b1 + b2;
    }
    typedef body<long double, long double, long double> bod;
    typedef body<long double, long double, long double, 0> bod_nf;
    typedef body<long double, long double, long double, 10> bod_10f;
    typedef body<long double, long double, long double, 100> bod_100f;
    typedef body<long double, long double, long double, 1000> bod_1000f;
}
#endif
