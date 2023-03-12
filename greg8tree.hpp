#ifndef GREG8TREE_H
#define GREG8TREE_H

#define RUNNING_MR_SUM

#ifndef __cplusplus
#error "greg8tree.hpp is a C++ header file only.\n"
#endif

#include "gregbod.hpp"

namespace gtd {
    class barnes_hut_error : public nbody_error {
    public:
        using nbody_error::nbody_error;
        barnes_hut_error() :
        nbody_error{"An error occurred with the creation and/or usage of a gtd::bh_cube<> object.\n"} {}
    };
    class invalid_body_placement : public nbody_error {
        char msg[96]{};
    public:
        invalid_body_placement() : msg{"Error: body does not fall within the current Barnes-Hut cube.\n"} {}
        invalid_body_placement(const char *message) {
            strcpy_c(this->msg, message);
        }
        invalid_body_placement(uint64_t id) {
            snprintf(msg, 96, "Error: body with ID = %" PRIu64" is not contained within the current Barnes-Hut cube.\n",
                     id);
        }
        const char *what() const noexcept override {
            return msg;
        }
    };
    class indistinguishable_bodies : public nbody_error {
        char msg[176]{};
    public:
        indistinguishable_bodies() : msg{"Error: body positions cannot be distinguished from one another.\n"} {}
        indistinguishable_bodies(const char *message) {
            strcpy_c(this->msg, message);
        }
        indistinguishable_bodies(uint64_t id1, uint64_t id2) {
            snprintf(msg, 176, "Error: body with ID = %" PRIu64" and body with ID = %" PRIu64" have positions that are "
                               "either equal or so close together that they cannot be told apart.\n", id1, id2);
        }
        const char *what() const noexcept override {
            return msg;
        }
    };
    template <isNumWrapper M = long double, isNumWrapper R = long double,
              isNumWrapper T = long double, uint64_t rF = 0>
    class bh_cube {
        using vec = vector3D<T>;
        using bod_t = body<M, R, T, rF>;
        using cube_t = bh_cube<M, R, T, rF>;
        vec _btm{}; // corner of cube with smallest x, y & z coordinates
        vec _top{}; // corner of cube with largest x, y & z coordinates
        vec _com{}; // centre of mass of all bodies contained within cube
        M _mass{}; // total mass of all bodies contained within cube
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
        // using mr_t = decltype(std::declval<M>()*std::declval<T>());
        vec mr_sum{};
#endif
#endif
        /* Only a leaf node (or cube) will have a _bod that is NOT a nullptr. */
        const bod_t *_bod; // pointer to body stored - only stores if leaf node
        cube_t *parent = nullptr; // parent cube
        /* Perhaps unconventionally, I am indexing the octants as follows: 0 -> (-,-,-), 1 -> (+,-,-),
         * 2 -> (-,+,-), 3 -> (+,+,-), 4 -> (-,-,+), 5 -> (+,-,+), 6 -> (-,+,+), 7 -> (+,+,+), */
        cube_t *_sub[8]{}; // 8 pointers to child cubes
        std::pair<vec, vec> corners[8]; // bottom and top corners of child cubes
        // cube_t *c1 = nullptr; // child cube 1
        // cube_t *c2 = nullptr; // etc...
        // cube_t *c3 = nullptr;
        // cube_t *c4 = nullptr;
        // cube_t *c5 = nullptr;
        // cube_t *c6 = nullptr;
        // cube_t *c7 = nullptr;
        // cube_t *c8 = nullptr;
        void set_corners() {
            /* Method to set the bounds of the 8 child cubes (even if they are never created). */
            T half_x = (_btm.x + _top.x)/2.0l;
            T half_y = (_btm.y + _top.y)/2.0l;
            T half_z = (_btm.z + _top.z)/2.0l;
            std::pair<vec, vec> *ptr = corners;
            ptr->first = _btm; // octant 0 bottom
            ptr->second.x = half_x; // octant 0 top
            ptr->second.y = half_y;
            ptr++->second.z = half_z;
            ptr->first.x = half_x; // octant 1 bottom
            ptr->first.y = _btm.y;
            ptr->first.z = _btm.z;
            ptr->second.x = _top.x; // octant 1 top
            ptr->second.y = half_y;
            ptr++->second.z = half_z;
            ptr->first.x = _btm.x; // octant 2 bottom
            ptr->first.y = half_y;
            ptr->first.z = _btm.z;
            ptr->second.x = half_x; // octant 2 top
            ptr->second.y = _top.y;
            ptr++->second.z = half_z;
            ptr->first.x = half_x; // octant 3 bottom
            ptr->first.y = half_y;
            ptr->first.z = _btm.z;
            ptr->second.x = _top.x; // octant 3 top
            ptr->second.y = _top.y;
            ptr++->second.z = half_z;
            ptr->first.x = _btm.x; // octant 4 bottom
            ptr->first.y = _btm.y;
            ptr->first.z = half_z;
            ptr->second.x = half_x; // octant 4 top
            ptr->second.y = half_y;
            ptr++->second.z = _top.z;
            ptr->first.x = half_x; // octant 5 bottom
            ptr->first.y = _btm.y;
            ptr->first.z = half_z;
            ptr->second.x = _top.x; // octant 5 top
            ptr->second.y = half_y;
            ptr++->second.z = _top.z;
            ptr->first.x = _btm.x; // octant 6 bottom
            ptr->first.y = half_y;
            ptr->first.z = half_z;
            ptr->second.x = half_x; // octant 6 top
            ptr->second.y = _top.y;
            ptr++->second.z = _top.z;
            ptr->first.x = half_x; // octant 7 bottom
            ptr->first.y = half_y;
            ptr->first.z = half_z;
            ptr->second = _top; // octant 7 top
        }
        void check_first_body() {
            if (_bod == nullptr)
                throw std::invalid_argument{"Error: pointer to gtd::body<M, R, T, rF> cannot be nullptr.\n"};
            if (!_bod->curr_pos.between(_btm, _top))
                throw invalid_body_placement{_bod->id};
            this->_mass = this->_bod->mass_;
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
            this->mr_sum = this->_mass*this->_bod->curr_pos;
#endif
#else
            this->_com = this->_bod->curr_pos;
#endif
        }
        void dispatch(const bod_t *_body) {
            std::pair<vec, vec> *ptr;
            cube_t *cptr = this;
            unsigned char i;
            // for (unsigned char i = 0; i < 8; ++i, ++ptr) {
#ifdef RUNNING_COM
            M prev_mass;
#endif
            while (true) {
                ptr = cptr->corners;
                for (i = 0; i < 8;) {
                    if (_body->curr_pos.between(ptr->first, ptr->second)) {
                        // if (_body->id == 0 && i == 7)
                        //     std::cout << "first: " << ptr->first << ", second: " << ptr->second << std::endl;
                        // if (_sub[i] == nullptr) {
                        //     _sub[i] = new bh_cube{this, ptr->first, ptr->second, _body};
                        // } else {
                        //     _sub[i]->add_body(_body);
                        // }
                        // return true;
                        // cptr->tot_mass += body->mass_;
                        if (cptr->_sub[i] == nullptr) {
                            cptr->_sub[i] = new bh_cube{cptr, ptr->first, ptr->second, _body};
                            // std::cout << "Created: " << cptr->_sub[i] << std::endl;
                            if (cptr->_bod != nullptr) { // means *cptr was a leaf node
                                if (_body->curr_pos == cptr->_bod->curr_pos)
                                    throw indistinguishable_bodies{_body->id, cptr->_bod->id};
                                _body = cptr->_bod;
                                cptr->_bod = nullptr;
                                ptr -= i;
                                i = 0;
                                continue;
                            }
                            return;
                        }
                        else {
                            cptr = cptr->_sub[i];
#ifdef RUNNING_COM
                            prev_mass = cptr->_mass;
#endif
                            cptr->_mass += _body->mass_;
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
                            cptr->mr_sum += _body->mass_*_body->curr_pos;
#endif
#else
                            cptr->_com = (cptr->_com*prev_mass + _body->mass_*_body->curr_pos)/cptr->_mass;
#endif
                            break;
                        }
                    }
                    // std::cout << "i: " << +i << std::endl;
                    ++i;
                    ++ptr;
                }
                if (i == 8) [[unlikely]] {
                    // std::cout << "bottom: " << cptr->_btm << ", top: " << cptr->_top << std::endl;
                    throw invalid_body_placement{_body->id}; // in case of f.p. error
                }
            }
            // throw invalid_body_placement{_body->id}; // in case of f.p. error
        }
#ifndef RUNNING_COM
    public:
        void calc_com() noexcept {
            if (this->_bod != nullptr) {
                this->_com = this->_bod->curr_pos;
                return;
            }
            this->_com = this->mr_sum/this->_mass; // calculate com of this cube outside loop
            cube_t *cptr = this;
            cube_t **sptr = this->_sub;
            unsigned char i = 0;
#ifdef RUNNING_MR_SUM
            while (true) {
                // cptr->_com = cptr->mr_sum/cptr->_mass;
                // sptr = cptr->_sub;
                while (i < 8) {
                    if (*sptr != nullptr) {
                        if ((*sptr)->_bod == nullptr) {
                            cptr = *sptr;
                            sptr = cptr->_sub;
                            i = 0;
                            continue;
                        } else {
                            (*sptr)->_com = (*sptr)->mr_sum/(*sptr)->_mass;
                        }
                    }
                    ++i;
                    ++sptr;
                }
                if (cptr->parent == nullptr)
                    break;
                cptr->_com = cptr->mr_sum/cptr->_mass;
                i = 1;
                sptr = cptr->parent->_sub;
                while (*sptr++ != cptr) ++i;
                cptr = cptr->parent;
            }
#else
#endif
        }
#endif
    public:
        // bh_cube() = default;
        // bh_cube(cube_t *parent_cube, const vec &btm, const vec &top) : parent{parent_cube},
        // _btm{btm}, _top{top} {this->set_corners();}
        // bh_cube(cube_t *parent_cube, vec &&btm, vec &&top) : parent{parent_cube},
        // _btm{std::move(btm)}, _top{std::move(_top)} {this->set_corners();}
        bh_cube(cube_t *parent_cube, const vec &btm, const vec &top, const bod_t *_body) : parent{parent_cube},
        _btm{btm}, _top{top}, _bod{_body} {
            this->set_corners();
            this->check_first_body();
        }
        bh_cube(cube_t *parent_cube, vec &&btm, vec &&top, const bod_t *_body) : parent{parent_cube},
        _btm{std::move(btm)}, _top{std::move(top)}, _bod{_body} {
            this->set_corners();
            this->check_first_body();
        }
        const bod_t *get_body() const noexcept {
            return this->_bod;
        }
        bool add_body(const bod_t *_body) {
            if (_body == nullptr) [[unlikely]]
                return false;
            if (!_body->curr_pos.between(this->_btm, this->_top)) [[unlikely]]
                throw invalid_body_placement{this->_bod->id};
            // if (this->_bod != nullptr) { // case for first split of cube into sub-cubes
            //     // if (!this->dispatch(this->_bod)) {
            //     //     throw invalid_body_placement{this->_bod->id};
            //     // }
            //     this->dispatch(this->_bod);
            //     this->_bod = nullptr;
            // }
            // this->_mass += _body->mass_;
            // if (!this->dispatch(_body))
            //     throw invalid_body_placement{_body->id};
#ifdef RUNNING_COM
            M prev_mass = this->_mass;
#endif
            this->_mass += _body->mass_;
#ifndef RUNNING_COM
#ifdef RUNNING_MR_SUM
            this->mr_sum += _body->mass_*_body->curr_pos;
#endif
#else
            this->_com = (this->_com*prev_mass + _body->mass_*_body->curr_pos)/this->_mass;
#endif
            this->dispatch(_body);
            return true;
        }
        vec cube_com() const noexcept {
            // this->calc_com();
            return this->_com;
        }
        ~bh_cube() {
            // std::cout << "Which dtor? this = " << this << std::endl;
            if (this->parent == nullptr) [[unlikely]] { // only root node should do the destructing
                cube_t *cptr = this;
                cube_t **sptr;
                unsigned char i = 0;
                // unsigned char count = 0;
                while (true) {
                    sptr = cptr->_sub;
                    loop_start:
                    for (; i < 8; ++i, ++sptr) {
                        if (*sptr != nullptr) {
                            if ((*sptr)->_bod != nullptr) { // case for leaf node
                                delete *sptr;
                                *sptr = nullptr;
                            }
                            else {
                                cptr = *sptr;
                                goto loop_end;
                            }
                        }
                    }
                    if (cptr->parent == nullptr) // got back up to the root node/cube
                        break;
                    if (i == 8) {
                        i = 0;
                        sptr = cptr->parent->_sub;
                        while (*sptr++ != cptr) ++i; // get *sptr to point to the next cube that needs processing
                        // i = cptr - *(cptr->parent->_sub);
                        cptr = cptr->parent;
                        delete cptr->_sub[i];
                        cptr->_sub[i++] = nullptr;
                        goto loop_start;
                    }
                    loop_end:
                    i = 0;
                }
                // while (cptr->_bod != nullptr) {
                //     sptr = cptr->_sub;
                //     for (i = 0; i < 8; ++i, ++sptr) {
                //         if (*sptr != nullptr) {
                //             cptr = *sptr;
                //             break;
                //         }
                //     }
                // }
            }
        }
        const cube_t &operator[](unsigned char index) {
            if (index > 7)
                throw std::invalid_argument{"gtd::bh_cube<M, R, T, rF>::operator[] error: index cannot be greater than "
                                            "7.\n"};
            if (this->_sub[index] == nullptr)
                throw std::invalid_argument{"Error: no gtd::bh_cube<> at index specified.\n"};
            return *(this->_sub[index]);
        }
        template <isNumWrapper m, isNumWrapper r, isNumWrapper t, uint64_t recFreq>
        friend std::ostream &operator<<(std::ostream &os, const bh_cube<m, r, t, recFreq> &cube);
        template <isNumWrapper, isNumWrapper, isNumWrapper, uint64_t>
        friend class bh_tree;
    };
    template <isNumWrapper m, isNumWrapper r, isNumWrapper t, uint64_t recFreq>
    std::ostream &operator<<(std::ostream &os, const bh_cube<m, r, t, recFreq> &cube) {
        /* The ONLY function in which I have allowed myself to use recursion (since it is not an important function).
         * Beware of stack overflows! */
        const bh_cube<m, r, t, recFreq> *const *ptr = cube._sub;
        os << "[gtd::bh_cube@" << &cube << ":mass=" << cube._mass << ",com=" << cube._com << ",body=";
        if (cube._bod == nullptr)
            os << "NULL";
        else
            os << *cube._bod;
        os << ",cubes:\n";
        for (unsigned char i = 0; i < 8; ++i, ++ptr) {
            os << "cube_" << +i << '=';
            if (*ptr == nullptr)
                os << "NULL";
            else
                os << '\n' << **ptr;
            os << '\n';
        }
        return os << ']';
    }
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T, uint64_t rF>
    class bh_tree {
        using typename bh_cube<M, R, T, rF>::cube_t;
        cube_t *_root = nullptr;
    public:
        bh_tree() = default;
        ~bh_tree() {
            delete _root;
        }
    };
    template <uint64_t recFreq>
    using cube = bh_cube<long double, long double, long double, recFreq>;
    typedef bh_cube<long double, long double, long double, 0> cube_0f;
    typedef bh_cube<long double, long double, long double, 1> cube_1f;
    typedef bh_cube<long double, long double, long double, 10> cube_10f;
    typedef bh_cube<long double, long double, long double, 100> cube_100f;
    typedef bh_cube<long double, long double, long double, 1000> cube_1000f;
}
#endif
