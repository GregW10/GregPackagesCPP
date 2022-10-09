//
// Created by mario on 08/10/2022.
//

#ifndef GREGBOD_H
#define GREGBOD_H

#include "gregvec.h"
#include <set>

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
        inline void add_pos_vel() {
            positions.push_back(curr_pos);
            velocities.push_back(curr_vel);
        }
        inline void add_pos_vel_ke() {
            add_pos_vel();
            set_ke();
            energies.push_back(curr_ke);
        }
        inline void set_ke() {
            auto mag = curr_vel.magnitude(); // best to push var. onto stack rather than call func. twice
            curr_ke = 0.5*mass*mag*mag; // this avoids the call to pow(mag, 2)
        }
        inline void check_mass() { // I prefer to throw an exception, rather than simply taking no action, to make it
            if (mass < M{0}) {     // clear where a negative quantity has attempted to be set
                throw negative_mass_error();
            }
        }
        inline void check_radius() {
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
        M get_mass() {
            return mass;
        }
        R get_radius() {
            return radius;
        }
        vector3D<T> get_pos() {
            return curr_pos;
        }
        vector3D<T> get_vel() {
            return curr_vel;
        }
        K get_ke() {
            return curr_ke;
        }
        inline void set_mass(const M &new_mass) {
            mass = new_mass;
            check_mass();
        }
        inline void set_radius(const R &new_radius) {
            radius = new_radius;
            check_radius();
        }
        inline void set_mass(M &&new_mass) {
            set_mass(new_mass);
        }
        inline void set_radius(R &&new_radius) {
            set_radius(new_radius);
        }
        void update(const vector3D<T> &new_position, const vector3D<T> &new_velocity) {
            curr_pos = new_position; // no checks to be performed, these new vecs can be anything
            curr_vel = new_velocity;
            add_pos_vel_ke();
        }
        void shift(const vector3D<T> &position_shift, const vector3D<T> &velocity_shift) {
            curr_pos += position_shift;
            curr_vel += velocity_shift;
            add_pos_vel_ke();
        }
        void update(vector3D<T> &&new_position, vector3D<T> &&new_velocity) {
            update(new_position, new_velocity);
        }
        void shift(vector3D<T> &&position_shift, vector3D<T> &&velocity_shift) {
            shift(position_shift, velocity_shift);
        }
        void clear() {
            positions.clear();
            velocities.clear();
            energies.clear();
            positions.push_back(curr_pos);
            velocities.push_back(curr_vel);
            energies.push_back(curr_ke);
        }
        std::vector<vector3D<T>> get_positions() {
            return {positions};
        }
        std::vector<vector3D<T>> get_velocity() {
            return {velocities};
        }
        std::vector<K> get_kinetic_energies() {
            return {energies};
        }
        auto begin() {
            return zip_ref(positions, velocities, energies).begin();
        }
        auto end() {
            return zip_ref(positions, velocities, energies).end();
        }
        auto cbegin() const {
            return zip_cref(positions, velocities, energies).cbegin();
        }
        auto cend() const {
            return zip_cref(positions, velocities, energies).cend();
        }
        auto rbegin() {
            return zip_ref(positions, velocities, energies).rbegin();
        }
        auto rend() {
            return zip_ref(positions, velocities, energies).rend();
        }
        auto crbegin() const {
            return zip_cref(positions, velocities, energies).crbegin();
        }
        auto crend() const {
            return zip_cref(positions, velocities, energies).crend();
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend std::ostream &operator<<(std::ostream &os, const body<m, r, t> &bod);
    };
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
    std::ostream &operator<<(std::ostream &os, const body<m, r, t> &bod) {
        return os << "[gtd::body@" << &bod << ":id=" << bod.id << ",m=" << bod.mass << ",r=" << bod.radius
                  << ",current_pos=(" << bod.curr_pos << "),current_vel=(" << bod.curr_vel << "),current_ke="
                  << bod.curr_ke << "]";
    }
    std::set<unsigned long long> body_counter::ids;
    unsigned long long body_counter::count = 0;
}
#endif
