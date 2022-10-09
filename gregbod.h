//
// Created by mario on 08/10/2022.
//

#ifndef GREGBOD_H
#define GREGBOD_H

#include "gregvec.h"
#include <set>
#include <fstream>
#include <map>

namespace gtd {
    class overlapping_bodies_error : public std::exception {
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
    class negative_mass_error : public std::exception {
        const char *message;
    public:
        negative_mass_error() : message("A gtd::body cannot have a negative mass.\n") {}
        explicit negative_mass_error(const char *msg) : message(msg) {}
        const char *what() const noexcept override {
            return message;
        }
    };
    class negative_radius_error : public std::exception {
        const char *message;
    public:
        negative_radius_error() : message("A gtd::body cannot have a negative radius.\n") {}
        explicit negative_radius_error(const char *msg) : message(msg) {}
        const char *what() const noexcept override {
            return message;
        }
    };
    /* body_counter was created because each instantiation of the body subclass with different template parameters is
     * actually a different class, so the count of bodies for each different class template instantiation would be
     * different */
    class body_counter {
        static unsigned long long count;
        static std::set<unsigned long long> ids;
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
    };
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
    class body : public body_counter {
        M mass;
        R radius;
        vector3D<T> curr_pos; // current position
        vector3D<T> curr_vel; // current velocity - have not included acc. as this should always be det. externally
        using K = decltype(0.5*mass*curr_vel.magnitude()*curr_vel.magnitude());
        K curr_ke; // current kinetic energy
        std::vector<vector3D<T>> positions; // these std::vectors will hold all the positions and velocities of the
        std::vector<vector3D<T>> velocities; // ... body as it moves
        std::vector<K> energies;
        mutable std::vector<std::tuple<const vector3D<T>&, const vector3D<T>&, const K&>> cref;
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
            curr_ke = 0.5*mass*mag*mag; // this avoids the call to pow(mag, 2)
        } // luckily, ^^^ 0.5 can be represented exaclty in binary (although mag probably won't be exact)
        void check_mass() { // I prefer to throw an exception, rather than simply taking no action, to make it clear
            if (mass < M{0}) { // ... where a negative quantity has attempted to be set
                throw negative_mass_error();
            }
        }
        void check_radius() {
            if (radius < R{0}) {
                throw negative_radius_error();
            }
        }
    public:
        body() : mass{1}, radius{1}, curr_ke{0}, curr_pos{}, curr_vel{} {add_pos_vel_ke(); check_mass();
            check_radius();}
        body(M &&body_mass, R &&body_radius) : mass{body_mass}, radius{body_radius}, curr_ke{0}, curr_pos{},
                                               curr_vel{} {add_pos_vel(); check_mass(); check_radius();}
        body(const M &body_mass, const R &body_radius) : mass{body_mass}, radius{body_radius}, curr_ke{0}, curr_pos{},
                                                         curr_vel{} {add_pos_vel(); check_mass(); check_radius();}
        body(M &&body_mass, R &&body_radius, vector3D<T> &&pos, vector3D<T> &&vel) :
        mass{body_mass}, radius{body_radius}, curr_pos{pos}, curr_vel{vel} {add_pos_vel_ke(); check_mass();
            check_radius();}
        body(const M &body_mass, const R &body_radius, const vector3D<T> &pos, const vector3D<T> &vel) :
        mass{body_mass}, radius{body_radius}, curr_pos{}, curr_vel{} {add_pos_vel_ke(); check_mass(); check_radius();}
        const M &get_mass() const noexcept {
            return mass;
        }
        const R &get_radius() const noexcept {
            return radius;
        }
        const vector3D<T> &get_pos() const noexcept {
            return curr_pos;
        }
        const vector3D<T> &get_vel() const noexcept {
            return curr_vel;
        }
        const K &get_ke() const noexcept {
            return curr_ke;
        }
        const vector3D<T> &pos_at(unsigned long long index) const {
            if (index >= positions.size()) {
                throw std::out_of_range("Requested position does not exist (index out of range).\n");
            }
            return positions[index];
        }
        const vector3D<T> &vel_at(unsigned long long index) const {
            if (index >= velocities.size()) {
                throw std::out_of_range("Requested velocity does not exist (index out of range).\n");
            }
            return velocities[index];
        }
        const K &ke_at(unsigned long long index) const {
            if (index >= energies.size()) {
                throw std::out_of_range("Requested kinetic energy does not exist (index out of range).\n");
            }
            return energies[index];
        }
        void set_mass(const M &new_mass) {
            mass = new_mass;
            check_mass();
        }
        void set_radius(const R &new_radius) {
            radius = new_radius;
            check_radius();
        }
        void set_mass(M &&new_mass) {
            set_mass(new_mass);
        }
        void set_radius(R &&new_radius) {
            set_radius(new_radius);
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
        vector3D<T> operator[](unsigned long long index) const { // returns a copy
            if (index >= positions.size()) {
                throw std::out_of_range("The specified position index is out of range.\n");
            }
            return positions[index];
        }
        bool write_trajectory_to_txt(String path = get_home_path<String>() + file_sep() + "Body_Trajectory_At_" +
                get_date_and_time() + ".txt", bool truncate = false, bool full_csv_style = false) const {
            if (!path.contains('.'))
                path.append_back(".txt");
            std::ofstream out(path.c_str(), truncate ? std::ios_base::trunc : std::ios_base::app);
            if (!out.good()) {
                return false;
            }
            if (full_csv_style) {
                out << "body_id,mass,radius\r\n" << body_counter::id << ',' << mass << ',' << radius << "\r\n"
                    << "position_x,position_y,position_z,velocity_x,velocity_y,velocity_z,kinetic_energy\r\n";
                for (const auto &[pos, vel, ke] : *this) {
                    out << pos.x << "," << pos.y << "," << pos.z << "," << vel.x << "," << vel.y << "," << vel.z << ","
                        << ke << "\r\n";
                }
                out << "\r\n";
                out.close();
                return true;
            }
            long long before = out.tellp();
            out << "Body ID: " << body_counter::id << ", Mass = " << mass << ", Radius = " << radius << "\n";
            long long after = out.tellp() - before;
            --after;
            for (long i = 0; i < after; ++i) {
                out.put('-');
            }
            out.put('\n');
            out << "position,velocity,kinetic_energy\n";
            for (const auto &[pos, vel, ke] : *this) {
                out << pos << "," << vel << "," << ke << "\n";
            }
            out << '\n';
            out.close();
            return true;
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend std::ostream &operator<<(std::ostream &os, const body<m, r, t> &bod);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &&b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &&b1, const body<m2, r2, t2> &b2);
        template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
        friend auto operator+(const body<m1, r1, t1> &&b1, const body<m2, r2, t2> &&b2);
    };
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
    std::ostream &operator<<(std::ostream &os, const body<m, r, t> &bod) {
        return os << "[gtd::body@" << &bod << ":id=" << bod.id << ",m=" << bod.mass << ",r=" << bod.radius
                  << ",current_pos=(" << bod.curr_pos << "),current_vel=(" << bod.curr_vel << "),current_ke="
                  << bod.curr_ke << "]";
    }
    template <isNumWrapper m1, isNumWrapper r1, isNumWrapper t1, isNumWrapper m2, isNumWrapper r2, isNumWrapper t2>
    auto operator+(const body<m1, r1, t1> &b1, const body<m2, r2, t2> &b2) {
        return body<decltype(b1.mass + b2.mass), long double, decltype(std::declval<t1>() + std::declval<t2>())>(b1.mass
        + b2.mass, cbrtl(b1.radius*b1.radius*b1.radius + b2.radius*b2.radius*b2.radius), b1.curr_pos + b2.curr_pos,
        b1.curr_vel + b2.curr_vel); // the returned body has a reset trajectory (since it is new!)
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
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
    class system {
        constexpr long double G = 0.0000000000667430; // Newtonian constant of Gravitation (m^3 kg^-1 s^-2)
        static const int two_body = 0; // static constants to determine which integration method to perform
        static const int euler = 1;
        static const int modified_euler = 2;
        static const int rk4 = 4;
        std::vector<body<M, R, T>> bodies;
        long double timestep = 1;
        unsigned long long iterations = 10000;
        using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        std::map<unsigned long long, std::vector<P>> pe; // map to store potential energies of bodies ({id, pe})
    };
    std::set<unsigned long long> body_counter::ids;
    unsigned long long body_counter::count = 0;
}
#endif
