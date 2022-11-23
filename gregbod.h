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
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
    class system;
    /* body_counter was created because each instantiation of the body subclass with different template parameters is
     * actually a different class, so the count of bodies for each different class template instantiation would be
     * different */
    class body_counter {
        static inline unsigned long long count = 0;
        static inline std::set<unsigned long long> ids;
        void set_id() {
            unsigned long long prev = 0;
            for (const auto &i : ids) {
                if (i - prev > 1) {
                    id = prev + 1;
                    return;
                }
            }
            id = *(ids.end()--) + 1;
        }
    protected:
        unsigned long long id; // is immutable after the object has been constructed
    public:
        body_counter() {
            id = count;
            if (ids.contains(id)) {
                set_id();
            }
            ids.insert(id);
            ++count;
        }
        body_counter(body_counter &&other) noexcept {
            this->id = other.id;
            other.id = -1; // would never, ever be reached
            while (ids.contains(other.id)) {
                --other.id;
            }
            ++count; // this will cause count to actually be the same as before the call to the constructor, since
        } // ... "other" is destroyed at the end of the call, reducing count by 1, so count must be incremented.
        body_counter(const body_counter &other) : body_counter() {} // as if constructing an empty object, since all
        unsigned long long get_id() const noexcept {                // ... ids must be unique
            return id;
        }
        static unsigned long long get_tot_bc() noexcept {
            return count;
        }
        body_counter &operator=(const body_counter &other) { // again, takes no effect as id must be unique
            return *this;
        }
        static void print() {
            for (const unsigned long long &i : ids) {
                printf("Id: %llu\n", i);
            }
        }
        ~body_counter() {
            --count;
            ids.erase(id);
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        friend class system;
    };
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double>
    class body : public body_counter {
    public:
        static inline bool allow_self_addition = false;
    protected:
        R radius;
        vector3D<T> curr_pos; // current position
    private:
        M mass_;
        vector3D<T> curr_vel; // current velocity - have not included acc. as this should always be det. externally
        vector3D<T> acc; // current acceleration of the body - only used within system class
        T pe; // current potential energy of the body - only used in system class since requires other bodies to calc.
        using K = decltype(0.5*mass_*curr_vel.magnitude()*curr_vel.magnitude());
        K curr_ke; // current kinetic energy
        /* Only quantities that represent the state of a body on its own are stored (below) - position, velocity, and
         * kinetic energy - whilst quantities that depend on the states of other bodies are not stored - acceleration
         * and potential energy. */
        std::vector<vector3D<T>> positions; // these std::vectors will hold all the positions and velocities of the
        std::vector<vector3D<T>> velocities; // ... body as it moves
        std::vector<K> energies;
        std::vector<T> r_terms; // these 3 std::vectors are used in system for numerous integration schemes
        std::vector<T> r_grad_terms;
        std::vector<T> v_grad_terms;
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
        void set_ke() {
            auto mag = curr_vel.magnitude(); // best to push var. onto stack rather than call func. twice
            curr_ke = 0.5*mass_*mag*mag; // this avoids the call to pow(mag, 2)
        } // luckily, ^^^ 0.5 can be represented exactly in binary (although mag probably won't be exact)
        void check_mass() { // I prefer to throw an exception, rather than simply taking no action, to make it clear
            if (mass_ < M{0}) { // ... where a negative quantity has attempted to be set
                throw negative_mass_error();
            }
        }
        void check_radius() {
            if (radius < R{0}) {
                throw negative_radius_error();
            }
        }
    public:
        body() : mass_{1}, radius{1}, curr_ke{0}, curr_pos{}, curr_vel{}
        {add_pos_vel_ke(); check_mass(); check_radius();}
        body(M &&body_mass, R &&body_radius) : mass_{std::move(body_mass)}, radius{std::move(body_radius)}, curr_ke{0},
        curr_pos{}, curr_vel{} {add_pos_vel(); check_mass(); check_radius();}
        body(const M &body_mass, const R &body_radius) : mass_{body_mass}, radius{body_radius}, curr_ke{0}, curr_pos{},
        curr_vel{} {add_pos_vel(); check_mass(); check_radius();}
        body(M &&body_mass, R &&body_radius, vector3D<T> &&pos, vector3D<T> &&vel) :
        mass_{std::move(body_mass)}, radius{std::move(body_radius)}, curr_pos{std::move(pos)}, curr_vel{std::move(vel)}
        {add_pos_vel_ke(); check_mass(); check_radius();}
        body(const M &body_mass, const R &body_radius, const vector3D<T> &pos, const vector3D<T> &vel) :
        mass_{body_mass}, radius{body_radius}, curr_pos{}, curr_vel{} {add_pos_vel_ke(); check_mass(); check_radius();}
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t>
        body(const body<m, r, t> &other) : mass_{other.mass_}, radius{other.radius}, curr_pos{other.curr_pos},
        curr_vel{other.curr_vel}, curr_ke{other.curr_ke}, acc{other.acc} {add_pos_vel();}
        body(const body<M, R, T> &other) : mass_{other.mass_}, radius{other.radius}, curr_pos{other.curr_pos},
        curr_vel{other.curr_vel}, curr_ke{other.curr_ke}, acc{other.acc} {add_pos_vel();}
        body(body<M, R, T> &&other) noexcept : body_counter{std::move(other)}, mass_{std::move(other.mass_)},
        radius{std::move(other.radius)}, curr_pos{std::move(other.curr_pos)}, curr_vel{std::move(other.curr_vel)},
        curr_ke{std::move(other.curr_ke)}, acc{std::move(other.acc)}, positions{std::move(other.positions)},
        velocities{std::move(other.velocities)}, energies{std::move(other.energies)}, cref{std::move(other.cref)} {}
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
        const vector3D<T> &acceleration() const noexcept {
            return acc;
        }
        const K &ke() const noexcept {
            return curr_ke;
        }
        const T &potential_energy() const noexcept {
            return pe;
        }
        const vector3D<T> &prev_pos_at(vec_s_t index) const {
            if (index >= positions.size()) {
                throw std::out_of_range("Requested position does not exist (index out of range).\n");
            }
            return positions[index];
        }
        const vector3D<T> &prev_vel_at(vec_s_t index) const {
            if (index >= velocities.size()) {
                throw std::out_of_range("Requested velocity does not exist (index out of range).\n");
            }
            return velocities[index];
        }
        const K &prev_ke_at(typename std::vector<K>::size_type index) const {
            if (index >= energies.size()) {
                throw std::out_of_range("Requested kinetic energy does not exist (index out of range).\n");
            }
            return energies[index];
        }
        void set_mass(const M &new_mass) {
            mass_ = new_mass;
            check_mass();
        }
        void set_radius(const R &new_radius) {
            radius = new_radius;
            check_radius();
        }
        void set_mass(const M &&new_mass) {
            set_mass(new_mass);
        }
        void set_radius(const R &&new_radius) {
            set_radius(new_radius);
        }
        void set_acc(const vector3D<T> &new_acc) {
            acc = new_acc;
        }
        void set_acc(vector3D<T> &&new_acc) noexcept {
            acc = std::move(new_acc);
        }
        void set_pe(const T &new_pe) {
            pe = new_pe;
        }
        void set_pe(T &&new_pe) noexcept {
            pe = std::move(new_pe);
        }
        void update(const vector3D<T> &new_position, const vector3D<T> &new_velocity) noexcept {
            curr_pos = new_position; // no checks to be performed, these new vecs can be anything
            curr_vel = new_velocity;
            add_pos_vel_ke();
        }
        void shift(const vector3D<T> &position_shift, const vector3D<T> &velocity_shift) noexcept {
            curr_pos += position_shift;
            curr_vel += velocity_shift;
            add_pos_vel_ke();
        }
        void update(vector3D<T> &&new_position, vector3D<T> &&new_velocity) noexcept {
            update(new_position, new_velocity);
        }
        void shift(vector3D<T> &&position_shift, vector3D<T> &&velocity_shift) noexcept {
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
        auto momentum() const noexcept {
            return mass_*curr_vel;
        }
        void reset_acc() noexcept {
            acc = T{0};
        }
        void reset(bool clear_trajectory = true) { // resets the body to initial setup and clears traj. if specified
            curr_pos = positions[0];
            curr_vel = velocities[0];
            curr_ke = energies[0];
            if (clear_trajectory) {
                clear();
            }
            add_pos_vel();
            energies.push_back(curr_ke);
        }
        void clear() {
            positions.clear();
            velocities.clear();
            energies.clear();
            cref.clear();
            add_pos_vel();
            energies.push_back(curr_ke); // no need to calculate ke again
        }
        std::vector<vector3D<T>> get_positions_cpy() const {
            return {positions};
        }
        std::vector<vector3D<T>> get_velocities_cpy() const {
            return {velocities};
        }
        std::vector<K> get_kinetic_energies_cpy() const {
            return {energies};
        }
        const std::vector<vector3D<T>> &get_positions() const {
            return positions;
        }
        const std::vector<vector3D<T>> &get_velocities() const {
            return velocities;
        }
        const std::vector<K> &get_kinetic_energies() const {
            return energies;
        }
        /* since one should not modify any past pos, vel or ke at all, and no current pos, vel or ke other than
         * through the update() or shift() methods, all the iteration methods below return const iterators */
        auto begin() const {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).cbegin();
            return cref.cbegin();
        }
        auto end() const {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).cend();
            return cref.cend();
        }
        auto cbegin() const {
            return begin();
        }
        auto cend() const {
            return end();
        }
        auto rbegin() const {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).crbegin();
            return cref.crbegin();
        }
        auto rend() const {
            if (cref.size() != positions.size())
                return (cref = zip_cref(positions, velocities, energies)).crend();
            return cref.crend();
        }
        auto crbegin() const {
            return rbegin();
        }
        auto crend() const {
            return rend();
        }
        bool trajectory_to_txt(std::ofstream &out, bool full_csv_style = true) const {
            if (!out.good())
                return false;
            unsigned long long count = 0;
            if (full_csv_style) {
                out << "body_id,mass,radius\r\n" << body_counter::id << ',' << mass_ << ',' << radius << "\r\n"
                    << "iteration,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy\r\n";
                for (const auto &[pos, vel, ke] : *this) {
                    out << count++ << ',' << pos.x << ',' << pos.y << ',' << pos.z << ',' << vel.x << ',' << vel.y
                        << ',' << vel.z << ',' << ke << "\r\n";
                }
                out << "\r\n";
                out.close();
                return true;
            }
            long long before = out.tellp();
            out << "Body ID: " << body_counter::id << ", Mass = " << mass_ << ", Radius = " << radius << "\n";
            long long after = out.tellp() - before;
            --after;
            for (long i = 0; i < after; ++i) {
                out.put('-');
            }
            out.put('\n');
            out << "iteration,position,velocity,kinetic_energy\n";
            for (const auto &[pos, vel, ke] : *this) {
                out << count++ << ',' << pos << ',' << vel << ',' << ke << '\n';
            }
            out << '\n';
            return true;
        }
        bool trajectory_to_txt(const String &path = def_path, bool truncate = false, bool full_csv_style = false) const{
            if (&path == &def_path)
                def_path.append_back(get_date_and_time()).append_back(".csv");
            std::ofstream out(path.c_str(), truncate ? std::ios_base::trunc : std::ios_base::app);
            bool ret = trajectory_to_txt(out, full_csv_style);
            out.close();
            if (&path == &def_path)
                def_path.erase_chars(def_path.get_length() - 29);
            return ret;
        }
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t>
        body<M, R, T> &operator=(const body<m, r, t> &other) {
            if (&other == this) {
                return *this;
            }
            this->mass_ = other.mass_;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
            this->acc = other.acc;
            this->pe = other.pe;
            clear();
            return *this;
        }
        body<M, R, T> &operator=(const body<M, R, T> &other) {
            if (&other == this) {
                return *this;
            }
            this->mass_ = other.mass_;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
            this->acc = other.acc;
            this->pe = other.pe;
            clear();
            return *this;
        }
        body<M, R, T> &operator=(body<M, R, T> &&other) noexcept {
            if (&other == this) {
                return *this;
            }
            this->mass_ = std::move(other.mass_);
            this->radius = std::move(other.radius);
            this->curr_pos = std::move(other.curr_pos);
            this->curr_vel = std::move(other.curr_vel);
            this ->curr_ke = std::move(other.curr_ke);
            this->acc = std::move(other.acc);
            this->pe = std::move(other.pe);
            this->positions = std::move(other.positions);
            this->velocities = std::move(other.velocities);
            this->energies = std::move(other.energies);
            this->cref = std::move(other.cref);
            return *this;
        }
        vector3D<T> operator[](vec_s_t index) const { // returns a copy
            if (index >= positions.size()) {
                throw std::out_of_range("The specified position index is out of range.\n");
            }
            return positions[index];
        }
        body<M, R, T> &operator+=(const body<M, R, T> &other) noexcept requires isConvertible<R, long double> {
            if (&other == this) {
                if (!allow_self_addition)
                    return *this;
                this->mass_ *= M{2};
                this->radius = cbrtl(2*this->radius*this->radius*this->radius);
                return *this; // clearly, for self-addition, position and velocity remain unchanged
            }
            this->curr_pos = com(*this, other);
            this->curr_vel = avg_vel(*this, other);
            this->mass_ += other.mass_;
            this->radius = cbrtl(this->radius*this->radius*this->radius + other.radius*other.radius*other.radius);
            this->add_pos_vel_ke();
            return *this;
        }
        virtual ~body() = default;
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend std::ostream &operator<<(std::ostream &os, const body<m, r, t> &bod);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend inline auto com(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend inline auto avg_vel(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &&b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &&b1, const body<m2, r2, t2> &b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &&b1, const body<m2, r2, t2> &&b2);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend system<m, r, t> operator+(const system<m, r, t> &sys1, const system<m, r, t> &sys2);
        template <isNumWrapper PosT, isNumWrapper DirT, isNumWrapper LenT>
        friend class ray;
        template <isNumWrapper M_U, isNumWrapper R_U, isNumWrapper T_U, isNumWrapper PosU, isNumWrapper DirU,
                isNumWrapper DistU, isNumWrapper LenU, isNumWrapper LumU>
        friend class astro_scene;
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend class system;
    };
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
    std::ostream &operator<<(std::ostream &os, const body<m, r, t> &bod) {
        return os << "[gtd::body@" << &bod << ":id=" << bod.id << ",m=" << +bod.mass_ << ",r=" << +bod.radius
                  << ",current_pos=(" << bod.curr_pos << "),current_vel=(" << bod.curr_vel << "),current_ke="
                  << bod.curr_ke << "]";
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    inline auto com(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2) { // centre of mass_ position
        return (b1.mass_*b1.curr_pos + b2.mass_*b2.curr_pos)/(b1.mass_ + b2.mass_);
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    inline auto avg_vel(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2) { /* average weighted velocity, taking
        conservation of momentum into account */
        return (b1.momentum() + b2.momentum())/(b1.mass_ + b2.mass_);
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2) { /* performs a "merging" of two bodies, with
        the new body having the sum of both masses, being at the centre-of-mass of both bodies, having a volume equal
        to the sum of both volumes (using the radii), and a new velocity such that momentum is conserved */
        /* self-addition is not considered here as the returned body is new (body has not been added onto itself) */
        return body<decltype(b1.mass_ + b2.mass_), long double, decltype(com(b1, b2))>(b1.mass_ + b2.mass_,
        cbrtl(b1.radius*b1.radius*b1.radius + b2.radius*b2.radius*b2.radius), com(b1, b2), avg_vel(b1, b2));
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &&b2) {
        return b1 + b2;
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    auto operator+(const body<m1, r1, t1> &&b1, const body<m2, r2, t2> &b2) {
        return b1 + b2;
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    auto operator+(const body<m1, r1, t1> &&b1, const body<m2, r2, t2> &&b2) {
        return b1 + b2;
    }
    template <isNumWrapper M = long double, isNumWrapper R = long double, isNumWrapper T = long double>
    class system { // G, below, has not been made static as it varies between objects, depending on units passed
    public:
        static constexpr int two_body = 1; // static constants to determine which integration method to perform
        static constexpr int euler = 2;
        static constexpr int modified_euler = 4;
        static constexpr int leapfrog_kdk = 8;
        static constexpr int leapfrog_dkd = 16;
        static constexpr int rk4 = 32;
        static constexpr int barnes_hut = 65536; // static constants to determine the force approx. method to use
        static constexpr int fast_multipole = 131072;
        static inline bool merge_if_overlapped = false;
    private:
        long double G = 66743; // Newtonian constant of Gravitation (* 10^15 m^3 kg^-1 s^-2)
        using vec_size_t = typename std::vector<body<M, R, T>>::size_type;
        using bod_t = body<M, R, T>;
        std::vector<bod_t> bods;
        long double dt = 1;
        unsigned long long iterations = 1000;
        long double prev_dt{}; // used in methods called after evolve(), since could be changed by setter
        unsigned long long prev_iterations{}; // same here
        // using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        std::map<unsigned long long, std::vector<T>> pe; // map to store potential energies of bodies ({id, pe})
        std::map<unsigned long long, std::vector<T>> energy; // total energy for each body at each iteration (KE + PE)
        std::vector<T> tot_pe; // total potential energy for the entire system at each iteration
        std::vector<T> tot_ke; // total kinetic energy for the entire system at each iteration
        std::vector<T> tot_e; // total energy for the entire system at each iteration (KE + PE)
        bool evolved = false;
        bool prog = false; // whether to show the progress of the evolution of the system
        static inline String def_path = get_home_path<String>() + FILE_SEP + "System_Trajectories_";
        static inline unsigned long long to_ull(String &&str) {
            if (!str.isnumeric()) {
                return -1;
            }
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
            G /= 1000000000000000; // correcting for the original G being *10^(15) times larger than it should be
        }
        void clear_bodies() { // makes sure the trajectories of all bodies are deleted
            for (auto &bod : bods) {
                bod.clear();
            }
        }
        void clear_bodies(vec_size_t &as_of) {
            vec_size_t size = bods.size();
            while (as_of < size) {
                bods[as_of++].clear();
            }
        }
        void clear_bodies(vec_size_t &&as_of) {
            clear_bodies(as_of);
        }
        void check_overlap() {
            vec_size_t size;
            const body<M, R, T> *outer;
            const body<M, R, T> *inner;
            vec_size_t i = 0;
            vec_size_t j;
            start:
            size = bods.size();
            for (; i < size; ++i) {
                outer = &bods[i];
                for (j = i + 1; j < size; ++j) {
                    inner = &bods[j];
                    if ((outer->curr_pos - inner->curr_pos).magnitude() < outer->radius + inner->radius) {
                        if (!merge_if_overlapped) {
                            String str = "The bodies with id=";
                            str.append_back(outer->id).append_back(" and id=").append_back(inner->id);
                            str.append_back(" that were added to this system object overlap.\n");
                            throw overlapping_bodies_error(str.c_str());
                        }
                        bods[i] += bods[j]; // merges the two overlapping bodies
                        bods.erase(bods.begin() + j); // thus the number of total bodies is reduced by 1
                        goto start;
                    }
                }
            }
        }
        void cumulative_acc(bod_t &b1, bod_t &b2) {
            vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
            long double &&r12_cubed_mag = (r12*r12*r12).magnitude();
            b1.acc += ((G*b2.mass_)/(r12_cubed_mag))*r12;
            b2.acc -= ((G*b1.mass_)/(r12_cubed_mag))*r12;
            auto &&pot_energy = -(G*b1.mass_*b2.mass_)/r12.magnitude();
            b1.pe += pot_energy;
            b2.pe += pot_energy;
        }
        // void cumulative_acc0(bod_t &b1, bod_t &b2) {
        //     vector3D<T> &&r12 = b2.curr_pos - b1.curr_pos;
        //     vector3D<T> &&r12_cubed = r12*r12*r12;
        //     b2.v_grad_terms[0] -= (b1.v_grad_terms[0] += ((G*b2.mass_)/(r12_cubed.magnitude()))*r12);
        //     b2.pe += (b1.pe -= (G*b1.mass_*b2.mass_)/r12.magnitude());
        // }
        //vector3D<T> kv2(const bod_t &bod) {
            //
        //}
        void take_euler_step() {
            for (bod_t &bod : bods) {
                bod.curr_pos = bod.curr_pos + bod.curr_vel*dt;
                bod.curr_vel = bod.curr_vel + bod.acc*dt;
                bod.add_pos_vel_ke();
            }
        }
        void take_modified_euler_step(const std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>>
                                      &predicted_vals) {
            unsigned long long count = 0;
            const long double &&step = dt/2;
            for (bod_t &bod : bods) {
                bod.curr_pos = bod.curr_pos + step*(bod.curr_vel + std::get<1>(predicted_vals[count]));
                bod.curr_vel = bod.curr_vel + step*(bod.acc + std::get<2>(predicted_vals[count++]));
                bod.add_pos_vel_ke();
            }
        }
        // void calc_initial_energy() {
        //     pe.clear();
        //     vec_size_t size = bods.size();
        //     for (bod_t &bod : bods)
        //         bod.pe = T{0};
        //     T total_pe = T{0};
        //     T total_ke = T{0};
        //     T total_e = T{0};
        //     for (unsigned long long outer = 0; outer < size; ++outer) {
        //         bod_t &ref = bods[outer];
        //         for (unsigned long long inner = outer + 1; inner < size; ++inner) {
        //             bods[inner].pe += (ref.pe -= (G*ref.mass_*bods[inner].mass_)/(ref.curr_pos -
        //                                                                           bods[inner].curr_pos).magnitude());
        //         }
        //         pe.emplace(ref.id, std::vector<T>{ref.pe});
        //         energy.emplace(ref.id, std::vector<T>{(T) (ref.curr_ke + ref.pe)});
        //         total_pe += ref.pe;
        //         total_ke += ref.curr_ke;
        //         total_e += ref.pe + ref.curr_ke;
        //     }
        //     tot_pe.push_back(total_pe);
        //     tot_ke.push_back(total_ke);
        //     tot_e.push_back(total_e);
        // }
        void create_energy_vectors() {
            vec_size_t size = bods.size();
            for (const bod_t &bod : bods) {
                pe.emplace(bod.id, std::vector<T>{});
                energy.emplace(bod.id, std::vector<T>{});
            }
        }
        void calc_final_energy() {
            // pe.clear();
            vec_size_t size = bods.size();
            for (bod_t &bod : bods)
                bod.pe = T{0};
            T total_pe{0};
            T total_ke{0};
            // T total_e{0};
            unsigned long long inner;
            for (unsigned long long outer = 0; outer < size; ++outer) {
                bod_t &ref = bods[outer];
                for (inner = outer + 1; inner < size; ++inner) {
                    bods[inner].pe += (ref.pe -= (G*ref.mass_*bods[inner].mass_)/(ref.curr_pos -
                                                                                  bods[inner].curr_pos).magnitude());
                }
                // pe.emplace(ref.id, std::vector<T>{ref.pe});
                // energy.emplace(ref.id, std::vector<T>{(T) (ref.curr_ke + ref.pe)});
                pe.find(ref.id)->second.push_back(ref.pe);
                energy.find(ref.id)->second.push_back(ref.curr_ke + ref.pe);
                total_pe += ref.pe;
                total_ke += ref.curr_ke;
                // total_e += ref.pe + ref.curr_ke;
            }
            total_pe /= T{2};
            tot_pe.push_back(total_pe);
            tot_ke.push_back(total_ke);
            tot_e.push_back(total_pe + total_ke);
        }
        static bool check_option(int option) noexcept {
            unsigned short loword = option & 0x0000ffff;
            unsigned short hiword = option >> 16;
            return (loword & (loword - 1)) == 0 || loword == 0 || loword > rk4 ||
                   (hiword & (hiword - 1)) == 0 || hiword > 2;
        }
    public:
        system(const std::initializer_list<bod_t> &list) : bods{list} {
            parse_units_format("M1:1,D1:1,T1:1");
            clear_bodies(0);
            check_overlap();
        }
        system(std::initializer_list<bod_t> &&list) : bods{std::move(list)} {
            parse_units_format("M1:1,D1:1,T1:1");
            clear_bodies(0);
            check_overlap();
        }
        explicit system(long double timestep = 1, unsigned long long num_iterations = 1000, bool show_progress = true,
                        const char *units_format = "M1:1,D1:1,T1:1") :
                        dt{timestep}, iterations{num_iterations}, prog{show_progress} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            parse_units_format(units_format);
        }
        system(const std::vector<body<M, R, T>> &bodies, long double timestep = 1,
               unsigned long long num_iterations = 1000, bool show_progress = true,
               const char *units_format = "M1:1,D1:1,T1:1") :
               bods{bodies}, dt{timestep}, iterations{num_iterations}, prog{show_progress} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            parse_units_format(units_format);
            clear_bodies(0);
            check_overlap();
        }
        system(std::vector<body<M, R, T>> &&bodies, long double timestep = 1,
               unsigned long long num_iterations = 1000, bool show_progress = true,
               const char *units_format = "M1:1,D1:1,T1:1") :
               bods{std::move(bodies)}, dt{timestep}, iterations{num_iterations}, prog{show_progress} {
            parse_units_format(units_format);
            clear_bodies(0);
            check_overlap();
        }
        system(const system<M, R, T> &other) :
        bods{other.bods}, dt{other.dt}, iterations{other.iterations}, prog{other.prog}, pe{other.pe}, G{other.G} {
            clear_bodies(0);
            check_overlap();
        }
        system(system<M, R, T> &&other) : bods{std::move(other.bods)}, dt{other.dt}, iterations{other.iterations},
        prog{other.prog}, pe{std::move(other.pe)}, G{other.G} {
            clear_bodies(0);
            check_overlap();
        }
        vec_size_t num_bodies() {
            return bods.size();
        }
        void set_progress(bool show_progress) noexcept {
            prog = show_progress;
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
            if (delta_t == 0)
                return false;
            dt = delta_t;
            return true;
        }
        bool set_timestep(const long double &&delta_t) noexcept {
            return set_timestep(delta_t);
        }
        system<M, R, T> &add_body(const body<M, R, T> &bod) {
            bods.push_back(bod);
            bods.back().clear(); // all bodies within a system object must start out without an evolution
            check_overlap();
            return *this;
        }
        system<M, R, T> &add_body(body<M, R, T> &&bod) {
            bods.push_back(std::move(bod));
            bods.back().clear();
            check_overlap();
            return *this;
        }
        system<M, R, T> &add_bodies(const std::vector<body<M, R, T>> &bodies) {
            if (bodies.size() == 0) {
                return *this;
            }
            vec_size_t index = bods.size();
            bods.insert(bods.end(), bodies.begin(), bodies.end());
            clear_bodies(index);
            check_overlap();
            return *this;
        }
        system<M, R, T> &add_bodies(std::vector<body<M, R, T>> &&bodies) {
            if (bodies.size() == 0) {
                return *this;
            }
            for (auto &b : bodies) {
                b.clear();
                bods.push_back(std::move(b));
            }
            check_overlap();
            return *this;
        }
        const body<M, R, T> &get_body(unsigned long long id) {
            for (const auto &b : *this) {
                if (b.id == id) {
                    return b;
                }
            }
            throw std::invalid_argument("The id passed does not correspond to any body present in this system "
                                        "object.\n");
        }
        bool remove_body(unsigned long long id) {
            auto end_it = bods.cend();
            for (typename std::vector<body<M, R, T>>::const_iterator it = bods.cbegin(); it < end_it; ++it) {
                if ((*it).id == id) {
                    bods.erase(it);
                    return true;
                }
            }
            return false;
        }
        void reset() {
            for (bod_t &bod : bods)
                bod.reset(true);
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
            for (bod_t &bod : bods)
                bod.clear();
            pe.clear();
            energy.clear();
            tot_pe.clear();
            tot_ke.clear();
            tot_e.clear();
            evolved = false;
        }
        bool evolve(int integration_method) {
            if (bods.size() == 0 || !check_option(integration_method))
                return false;
            time_t start = time(nullptr);
            unsigned long long steps = 0;
            if (evolved)
                this->clear_evolution();
            else
                this->create_energy_vectors();
            // this->calc_initial_energy();
            vec_size_t num_bods = bods.size();
            T total_pe;
            T total_ke;
            T total_e;
            if ((integration_method & two_body) == two_body) {
                if (bods.size() != 2)
                    throw two_body_error();
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                    return true;
                }
                return true;
            }
            else if ((integration_method & euler) == euler) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {
                    unsigned long long inner;
                    while (steps++ < iterations) {
                        total_pe = T{0};
                        total_ke = T{0};
                        // total_e = T{0};
                        for (bod_t &bod : bods) {
                            bod.acc = vector3D<T>::zero;
                            bod.pe = T{0};
                        }
                        for (unsigned long long outer = 0; outer < num_bods; ++outer) {
                            bod_t &ref = bods[outer];
                            for (inner = outer + 1; inner < num_bods; ++inner) {
                                cumulative_acc(ref, bods[inner]);
                            }
                            pe.find(ref.id)->second.push_back(ref.pe);
                            energy.find(ref.id)->second.push_back((T) (ref.curr_ke + ref.pe));
                            total_pe += ref.pe;
                            total_ke += ref.curr_ke;
                            // total_e += ref.pe + ref.curr_ke;
                        }
                        total_pe /= T{2}; // pairs of particles were counted twice, so must divide by 2
                        tot_pe.push_back(total_pe);
                        tot_ke.push_back(total_ke);
                        tot_e.push_back(total_pe + total_ke);
                        this->take_euler_step();
                        if (prog) // re-evaluating this within the loop is not great, but I don't want to dup. code
                            printf("Iteration %llu/%llu\r", steps, iterations);
                    }
                    if (prog) {
                        time_t total = time(nullptr) - start;
                        printf("\n--------------------Done--------------------\n"
                               "Euler method time elapsed: %zu second%c\n", total, "s"[total == 1]);
                    }
                }
            }
            else if ((integration_method & modified_euler) == modified_euler) {
                /* The below std::vector will hold the euler-predicted values for position, velocity and acceleration,
                 * respectively. These are subsequently used to determine the corrected values. */
                std::vector<std::tuple<vector3D<T>, vector3D<T>, vector3D<T>>> predicted{num_bods};
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                    return true;
                }
                else {
                    unsigned long long outer;
                    unsigned long long inner;
                    while (steps++ < iterations) {
                        total_pe = T{0};
                        total_ke = T{0};
                        // total_e = T{0};
                        for (bod_t &bod : bods) {
                            bod.acc = vector3D<T>::zero;
                            bod.pe = T{0};
                        }
                        for (outer = 0; outer < num_bods; ++outer) {
                            bod_t &ref = bods[outer];
                            for (inner = outer + 1; inner < num_bods; ++inner) {
                                cumulative_acc(ref, bods[inner]);
                            }
                            pe.find(ref.id)->second.push_back(ref.pe);
                            energy.find(ref.id)->second.push_back((T) (ref.curr_ke + ref.pe));
                            total_pe += ref.pe;
                            total_ke += ref.curr_ke;
                            // total_e += ref.pe + ref.curr_ke;
                        }
                        total_pe /= T{2};
                        tot_pe.push_back(total_pe);
                        tot_ke.push_back(total_ke);
                        tot_e.push_back(total_pe + total_ke);
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
                        if (prog)
                            printf("Iteration %llu/%llu\r", steps, iterations);
                    }
                    if (prog) {
                        time_t total = time(nullptr) - start;
                        printf("\n--------------------Done--------------------\n"
                               "Modified Euler method, time elapsed: %zu second%c\n", total, "s"[total == 1]);
                    }
                }
            }
            else if ((integration_method & leapfrog_kdk) == leapfrog_kdk) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                    return true;
                }
                return true;
            }
            else if ((integration_method & leapfrog_dkd) == leapfrog_dkd) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {

                }
            }
            else if ((integration_method & rk4) == rk4) {
                if ((integration_method & barnes_hut) == barnes_hut) {
                    // lots of stuff here
                }
                else {

                }
            }
            evolved = true;
            this->calc_final_energy();
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
            unsigned long long i;
            for (const bod_t &bod : bods) {
                out << "body_id,mass,radius\r\n" << bod.id << ',' << bod.mass_ << ',' << bod.radius << "\r\n"
                    << "time_elapsed,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy,"
                       "potential_energy,total_energy\r\n";
                for (i = 0; i <= prev_iterations; ++i) {
                    out << i*prev_dt << ',' << bod.positions[i].x << ',' << bod.positions[i].y << ','
                        << bod.positions[i].z << ',' << bod.velocities[i].x << ',' << bod.velocities[i].y << ','
                        << bod.velocities[i].z << ',' << bod.energies[i] << ',' << pe.find(bod.id)->second[i] << ','
                        << energy.find(bod.id)->second[i] << "\r\n";
                }
                out << ",,,,,,,,,\r\n";
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
            return true;
        }
        const body<M, R, T> &operator[](vec_size_t index) {
            if (index >= bods.size()) {
                throw std::out_of_range("The requested body does not exist (index out of range).\n");
            }
            return bods[index];
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend std::ostream &operator<<(std::ostream &os, const system<m, r, t> &sys);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend system<m, r, t> operator+(const system<m, r, t> &sys1, const system<m, r, t> &sys2);
    };
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
    std::ostream &operator<<(std::ostream &os, const system<M, R, T> &sys) {
        typename std::vector<body<M, R, T>>::size_type count = sys.bods.size();
        os << "[gtd::system@" << &sys << ",num_bodies=" << count;
        if (count == 0) {
            return os << ']';
        }
        os << ",bodies:";
        count = 0;
        for (const body<M, R, T> &b : sys) {
            os << "\n body_" << count++ << '=' << b;
        }
        os << ']';
        return os;
    }
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
    system<M, R, T> operator+(const system<M, R, T> &sys1, const system<M, R, T> &sys2) {
        typename std::vector<body<M, R, T>>::size_type size1 = sys1.bods.size();
        typename std::vector<body<M, R, T>>::size_type size2 = sys2.bods.size();
        if ((size1 == 0 && size2 == 0) || sys1.G != sys2.G) {
            return {{}};
        }
        if (size1 == 0) {
            return {sys2};
        }
        if (size2 == 0) {
            return {sys1};
        }
        system<M, R, T> ret_sys(sys1.bods, (sys1.dt + sys2.dt)/2.0l, (sys1.iterations + sys2.iterations)/2,
                                sys1.prog && sys2.prog);
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
    typedef body<long double, long double, long double> bod;
    typedef system<long double, long double, long double> sys;
}
#endif
