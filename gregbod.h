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
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
    class system;
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
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        friend class system;
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
        typedef typename std::vector<vector3D<T>>::size_type vec_s_t;
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
        body() : mass{1}, radius{1}, curr_ke{0}, curr_pos{}, curr_vel{} {add_pos_vel_ke(); check_mass();check_radius();}
        body(M &&body_mass, R &&body_radius) : mass{body_mass}, radius{body_radius}, curr_ke{0}, curr_pos{},
                                               curr_vel{} {add_pos_vel(); check_mass(); check_radius();}
        body(const M &body_mass, const R &body_radius) : mass{body_mass}, radius{body_radius}, curr_ke{0}, curr_pos{},
                                                         curr_vel{} {add_pos_vel(); check_mass(); check_radius();}
        body(M &&body_mass, R &&body_radius, vector3D<T> &&pos, vector3D<T> &&vel) :
        mass{body_mass}, radius{body_radius}, curr_pos{pos}, curr_vel{vel} {add_pos_vel_ke(); check_mass();
            check_radius();}
        body(const M &body_mass, const R &body_radius, const vector3D<T> &pos, const vector3D<T> &vel) :
        mass{body_mass}, radius{body_radius}, curr_pos{}, curr_vel{} {add_pos_vel_ke(); check_mass(); check_radius();}
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t>
        body(const body<m, r, t> &other) : mass{other.mass}, radius{other.radius}, curr_pos{other.curr_pos},
        curr_vel{other.curr_vel}, curr_ke{other.curr_ke} {add_pos_vel();}
        body(const body<M, R, T> &other) : mass{other.mass}, radius{other.radius}, curr_pos{other.curr_pos},
        curr_vel{other.curr_vel}, curr_ke{other.curr_ke} {add_pos_vel();}
        body(body<M, R, T> &&other) : body_counter{std::move(other)}, mass{other.mass}, radius{other.radius},
        curr_pos{other.curr_pos}, curr_vel{other.curr_vel}, curr_ke{other.curr_ke},
        positions{std::move(other.positions)}, velocities{std::move(other.velocities)},
        energies{std::move(other.energies)}, cref{std::move(other.cref)} {}
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
        const vector3D<T> &pos_at(vec_s_t index) const {
            if (index >= positions.size()) {
                throw std::out_of_range("Requested position does not exist (index out of range).\n");
            }
            return positions[index];
        }
        const vector3D<T> &vel_at(vec_s_t index) const {
            if (index >= velocities.size()) {
                throw std::out_of_range("Requested velocity does not exist (index out of range).\n");
            }
            return velocities[index];
        }
        const K &ke_at(typename std::vector<K>::size_type index) const {
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
        void set_mass(const M &&new_mass) {
            set_mass(new_mass);
        }
        void set_radius(const R &&new_radius) {
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
        template <isConvertible<M> m, isConvertible<R> r, isConvertible<T> t>
        body<M, R, T> &operator=(const body<m, r, t> &other) {
            if (&other == this) {
                return *this;
            }
            this->mass = other.mass;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
            clear();
            return *this;
        }
        body<M, R, T> &operator=(const body<M, R, T> &other) {
            if (&other == this) {
                return *this;
            }
            this->mass = other.mass;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
            clear();
            return *this;
        }
        body<M, R, T> &operator=(body<M, R, T> &&other) {
            if (&other == this) {
                return *this;
            }
            this->mass = other.mass;
            this->radius = other.radius;
            this->curr_pos = other.curr_pos;
            this->curr_vel = other.curr_vel;
            this ->curr_ke = other.curr_ke;
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
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend system<m, r, t> operator+(const system<m, r, t> &sys1, const system<m, r, t> &sys2);
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t>
        friend class system;
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
    class system { // G, below, has not been made static as it varies between objects, depending on units passed
        long double G = 66743; // Newtonian constant of Gravitation (* 10^15 m^3 kg^-1 s^-2)
        static constexpr int two_body = 0; // static constants to determine which integration method to perform
        static constexpr int euler = 1;
        static constexpr int modified_euler = 2;
        static constexpr int rk4 = 4;
        using vec_size_t = typename std::vector<body<M, R, T>>::size_type;
        std::vector<body<M, R, T>> bods;
        long double dt;
        unsigned long long iterations;
        using P = decltype((G*std::declval<M>()*std::declval<M>())/std::declval<T>()); // will be long double
        std::map<unsigned long long, std::vector<P>> pe; // map to store potential energies of bodies ({id, pe})
        bool prog; // whether to show the progress of the evolution of the system
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
        void check_overlap() {
            vec_size_t size = bods.size();
            const body<M, R, T> *outer;
            const body<M, R, T> *inner;
            for (vec_size_t i = 0; i < size; ++i) {
                outer = &bods[i];
                for (vec_size_t j = i + 1; j < size; ++j) {
                    inner = &bods[j];
                    if ((outer->curr_pos - inner->curr_pos).magnitude() < outer->radius + inner->radius) {
                        String str = "The bodies with id=";
                        str.append_back(outer->id).append_back(" and id=").append_back(inner->id);
                        str.append_back(" that were added to this system object overlap.\n");
                        throw overlapping_bodies_error(str.c_str());
                    }
                }
            }
        }
    public:
        system(const std::vector<body<M, R, T>> &bodies, long double timestep = 1,
               unsigned long long num_iterations = 1000, bool show_progress = true,
               const char *units_format = "M1:1,D1:1,T1:1") :
               bods{bodies}, dt{timestep}, iterations{num_iterations}, prog{show_progress} {
            /* units_format is a string with 3 ratios: it specifies the ratio of the units used for mass, distance and
             * time to kg, metres and seconds (SI units), respectively. This means any units can be used. */
            parse_units_format(units_format);
            clear_bodies();
            check_overlap();
        }
        system(std::vector<body<M, R, T>> &&bodies, long double timestep = 1,
               unsigned long long num_iterations = 1000, bool show_progress = true,
               const char *units_format = "M1:1,D1:1,T1:1") :
               bods{std::move(bodies)}, dt{timestep}, iterations{num_iterations}, prog{show_progress} {
            parse_units_format(units_format);
            clear_bodies();
            check_overlap();
        }
        system(const system<M, R, T> &other) :
        bods{other.bods}, dt{other.dt}, iterations{other.iterations}, prog{other.prog}, pe{other.pe}, G{other.G} {
            clear_bodies();
            check_overlap();
        }
        system(system<M, R, T> &&other) : bods{std::move(other.bods)}, dt{other.dt}, iterations{other.iterations},
        prog{other.prog}, pe{std::move(other.pe)}, G{other.G} {
            clear_bodies();
            check_overlap();
        }
        vec_size_t num_bodies() {
            return bods.size();
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
            return *this;
        }
        system<M, R, T> &add_bodies(std::vector<body<M, R, T>> &&bodies) {
            if (bodies.size() == 0) {
                return *this;
            }
            for (auto &b : bodies) {
                bods.push_back(std::move(b));
            }
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
        for (typename std::vector<body<M, R, T>>::size_type i = 0; i < size1; ++i) {
            for (typename std::vector<body<M, R, T>>::size_type j = 0; j < size2; ++j) {
                if ((sys1.bods[i].curr_pos-sys2.bods[j].curr_pos).magnitude() < sys1.bods[i].radius+sys2.bods[j].radius)
                    throw overlapping_bodies_error();
            }
        }
        system<M, R, T> ret_sys(sys1.bods, (sys1.dt + sys2.dt)/2.0l, (sys1.iterations + sys2.iterations)/2,
                                sys1.prog && sys2.prog);
        ret_sys.G = sys1.G;
        ret_sys.add_bodies(sys2.bods);
        return ret_sys;
    }
    std::set<unsigned long long> body_counter::ids;
    unsigned long long body_counter::count = 0;
}
#endif
