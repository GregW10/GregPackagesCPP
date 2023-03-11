#ifndef GREG8TREE_H
#define GREG8TREE_H

#define RUNNING_MR_SUM

#ifndef __cplusplus
#error "greg8tree.hpp is a C++ header file only.\n"
#endif

#include "gregbod.hpp"

namespace gtd {
    class invalid_body_placement : public nbody_error {
        char msg[96]{};
    public:
        invalid_body_placement() : msg{"Error: body does not fall within current Barnes-Hut cube.\n"} {}
        invalid_body_placement(const char *message) : msg{message} {}
        invalid_body_placement(uint64_t id) {
            snprintf(msg, 96, "Error: body with ID = %" PRIu64 " is not contained within current Barnes-Hut cube.\n",
                     id);
        }
        const char *what() const noexcept override {
            return msg;
        }
    };
    bool between(const vector3D<U> &v1, const vector3D<V> &v2) {
        /* Tests if the vector lies within the 3D cube spanned by v1 & v2 (v1 == bottom corner, v2 == top corner) */
        return v1.x <= this->x && this->x < v2.x &&
               v1.y <= this->y && this->y < v2.y &&
               v1.z <= this->z && this->z < v2.z &&;
    }
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
#ifdef RUNNING_MR_SUM
        // using mr_t = decltype(std::declval<M>()*std::declval<T>());
        vec mr_sum{};
#endif
        /* Only a leaf node (or cube) will have a _bod that is NOT a nullptr. */
        bod_t *_bod; // pointer to body stored - only stores if leaf node
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
            T half_x = (btm_corner.x + top_corner.x)/2.0l;
            T half_y = (btm_corner.y + top_corner.y)/2.0l;
            T half_z = (btm_corner.z + top_corner.z)/2.0l;
            std::pair<vec, vec> *ptr = corners;
            ptr->first = _btm;
            ptr->second.x = half_x;
            ptr->second.y = half_y;
            ptr++->second.z = half_z;
            ptr->first.x = half_x;
            ptr->first.y = _btm.y;
            ptr->first.z = _btm.z;
            ptr->second.x = _top.x;
            ptr->second.y = half_y;
            ptr++->second.z = half_z;
            ptr->first.x = _btm.x;
            ptr->first.y = half_y;
            ptr->first.z = _btm.z;
            ptr->second.x = half_x;
            ptr->second.y = _top.y;
            ptr++->second.z = half_z;
            ptr->first.x = half_x;
            ptr->first.y = half_y;
            ptr->first.z = _btm.z;
            ptr->second.x = _top.x;
            ptr->second.y = _top.y;
            ptr++->second.z = half_z;
            ptr->first.x = _btm.x;
            ptr->first.y = _btm.y;
            ptr->first.z = half_z;
            ptr->second.x = half_x;
            ptr->second.y = half_y;
            ptr++->second.z = _top.z;
            ptr->first.x = half_x;
            ptr->first.y = _btm.y;
            ptr->first.z = half_z;
            ptr->second.x = _top.x;
            ptr->second.y = half_y;
            ptr++->second.z = _top.z;
            ptr->first.x = _btm.x;
            ptr->first.y = half_y;
            ptr->first.z = half_z;
            ptr->second.x = half_x;
            ptr->second.y = _top.y;
            ptr++->second.z = _top.z;
            ptr->first.x = half_x;
            ptr->first.y = half_y;
            ptr->first.z = half_z;
            ptr->second = _top;
        }
        void check_first_body() {
            if (_bod == nullptr)
                throw std::invalid_argument{"Error: pointer to gtd::body<M, R, T, rF> cannot be nullptr.\n"};
            if (!_bod->between(_btm, _top))
                throw invalid_body_placement{_bod->id};
            this->_mass = this->_bod->mass_;
#ifdef RUNNING_MR_SUM
            this->mr_sum = this->_mass*this->_bod->curr_pos;
#endif
        }
        void dispatch(const bod_t *_body) {
            std::pair<vec, vec> *ptr;
            cube_t *cptr = this;
            unsigned char i;
            // for (unsigned char i = 0; i < 8; ++i, ++ptr) {
            while (true) {
                ptr = cptr->corners;
                for (i = 0; i < 8; ++i, ++ptr) {
                    if (_body->curr_pos.between(ptr->first, ptr->second)) {
                        // if (_sub[i] == nullptr) {
                        //     _sub[i] = new bh_cube{this, ptr->first, ptr->second, _body};
                        // } else {
                        //     _sub[i]->add_body(_body);
                        // }
                        // return true;
                        // cptr->tot_mass += body->mass_;
                        if (cptr->_sub[i] == nullptr) {
                            cptr->_sub[i] = new bh_cube{this, ptr->first, ptr->second, _body};
                            if (cptr->_bod != nullptr) { // means *cptr was a leaf node
                                _body = cptr->_bod;
                                cptr->_bod = nullptr;
                                i = 0;
                                continue;
                            }
                            return;
                        }
                        else {
                            cptr = cptr->_sub[i];
                            break;
                        }
                    }
                }
                cptr->_mass += body->mass_;
#ifdef RUNNING_MR_SUM
                this->mr_sum += body->mass_*body->curr_pos;
#endif
            }
            throw invalid_body_placement{_body->id}; // in case of f.p. error
        }
        void calc_com() {
            if (this->_bod != nullptr) {
                this->_com = this->_bod->curr_pos;
                return;
            }
            cube_t *cptr = this;
            cube_t **sptr;
            unsigned char i;
#ifdef RUNNING_MR_SUM
            while (true) {
                cptr->_com = cptr->mr_sum/cptr->_mass;
                sptr = cptr->_sub;
                for (i = 0; i < 8; ++i, ++sptr) {
                    if (*sptr != nullptr) {

                    }
                }
            }
#else
#endif
        }
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
            if (!_body->between(this->_btm, this->_top)) [[unlikely]]
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
            this->dispatch(_body);
            return true;
        }
        vec cube_com() {
            this->calc_com();
            return this->_com;
        }
        ~bh_cube() {
            if (this->parent == nullptr) [[unlikely]] { // only root node should do the destructing
                cube_t *cptr = this;
                cube_t **sptr;
                unsigned char i = 0;
                // unsigned char count = 0;
                while (true) {
                    sptr = cptr->_sub;
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
                        i = cptr - *(cptr->parent->_sub);
                        cptr = cptr->parent;
                        delete cptr->_sub[i];
                        cptr->_sub[i++] = nullptr;
                        continue;
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
    };
    class bh_tree {
        bh_cube *_root = nullptr;
    public:
        bh_tree() = default;
        ~bh_tree() {
            delete _root;
        }
    };
}
#endif
