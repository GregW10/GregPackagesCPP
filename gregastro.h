//
// Created by Gregor Hartl Watters on 04/11/2022.
//

#ifndef GREGASTRO_H
#define GREGASTRO_H

#include <functional>
#include <random>

#include "gregbod.h"
#include "gregbmp.h"

namespace gtd {
    class no_direction : public std::invalid_argument {
    public:
        no_direction() : std::invalid_argument{"Zero vector error.\n"} {}
        explicit no_direction(const char *message) : std::invalid_argument{message} {}
    };
    class zero_image_distance : public std::invalid_argument {
    public:
        zero_image_distance() : std::invalid_argument{"The image distance of a camera cannot be zero.\n"} {}
        explicit zero_image_distance(const char *message) : std::invalid_argument{message} {}
    };
    struct image_dimensions {
        unsigned int x;
        unsigned int y;
    };
    template <isNumWrapper PosT = long double, isNumWrapper DirT = long double, isNumWrapper LenT = long double>
    class ray : public vector3D<DirT> {
        /* A ray object represents a ray in simulated 3D space - whether this be an actual ray of light (from a light
         * source), or a "line-of-sight" ray, such as those shot out from a camera object (see below). A ray has an
         * origin (its position) a direction, and a parameter named 'l', which represents the length of the ray. A ray
         * with length zero is taken as one which continues on to infinity. Since the ray class is a child class of
         * vector3D<DirT>, the internal parent vector3D<DirT> object is made to be the direction, whilst a ray's
         * position is represented by a separate vector3D<PosT> object. */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        static inline const vector3D<DirT> zero{};
        vector3D<PosT> pos; // origin of the ray
        LenT l{0}; // length of the ray
        ray<PosT, DirT, LenT> &calc() {
            if (*this == zero)
                throw no_direction("A ray's direction cannot be a zero vector, "
                                   "or else the ray would have nowhere to point to.\n");
            this->normalise();
            return *this;
        }
        template <isNumWrapper T>
        static inline const T &&minimum(const T &&first, const T &&second) {
            return first <= second ? std::move(first) : std::move(second);
        }
    public:
        ray() : vector3D<DirT>{0, 0, -1} {calc();}
        ray(const vector2D<DirT> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        ray(vector2D<DirT> &&direction) noexcept : vector3D<DirT>(std::move(direction)) {calc();}
        ray(const vector3D<DirT> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        ray(vector3D<DirT> &&direction) noexcept : vector3D<DirT>(std::move(direction)) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &&direction) noexcept : vector3D<DirT>(direction) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &direction) noexcept : vector3D<DirT>(direction) {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &&direction) noexcept : vector3D<DirT>(direction) {calc();}
        ray(const DirT &x_component, const DirT &y_component, const DirT &z_component) noexcept :
                vector3D<DirT>{x_component, y_component, z_component} {calc();}
        ray(DirT &&x_component, DirT &&y_component, DirT &&z_component) noexcept :
        vector3D<DirT>{std::move(x_component), std::move(y_component), std::move(z_component)} {calc();}
        ray(const DirT &x_component, const DirT &y_component) noexcept :
        vector3D<DirT>{x_component, y_component} {calc();}
        ray(DirT &&x_component, DirT &&y_component) noexcept :
        vector3D<DirT>{std::move(x_component), std::move(y_component)} {calc();}
        ray(const vector2D<DirT> &origin, const vector2D<DirT> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        ray(vector2D<DirT> &&origin, vector2D<DirT> &&direction) noexcept :
        vector3D<DirT>(std::move(direction)), pos{std::move(origin)} {calc();}
        ray(const vector3D<DirT> &origin, const vector3D<DirT> &direction) noexcept :
        vector3D<DirT>(direction) {calc();}
        ray(vector3D<DirT> &&origin, vector3D<DirT> &&direction) noexcept :
        vector3D<DirT>(std::move(direction)), pos{std::move(origin)} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &origin, const vector2D<U> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector2D<U> &&origin, const vector2D<U> &&direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &origin, const vector3D<U> &direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        template <typename U> requires isConvertible<U, DirT>
        ray(const vector3D<U> &&origin, const vector3D<U> &&direction) noexcept :
        vector3D<DirT>(direction), pos{origin} {calc();}
        ray(const DirT &x_pos_comp, const DirT &y_pos_comp, const DirT &z_pos_comp,
            const DirT &x_dir_comp, const DirT &y_dir_comp, const DirT &z_dir_comp) noexcept :
                vector3D<DirT>{x_dir_comp, y_dir_comp, z_dir_comp}, pos{x_pos_comp, y_pos_comp, z_pos_comp} {calc();}
        ray(DirT &&x_pos_comp, DirT &&y_pos_comp, DirT &&z_pos_comp,
            DirT &&x_dir_comp, DirT &&y_dir_comp, DirT &&z_dir_comp) noexcept :
                vector3D<DirT>{std::move(x_dir_comp), std::move(y_dir_comp), std::move(z_dir_comp)},
                pos{std::move(x_pos_comp), std::move(y_pos_comp), std::move(z_pos_comp)} {calc();}
        ray(const DirT &x_pos_comp, const DirT &y_pos_comp, const DirT &x_dir_comp, const DirT &y_dir_comp) noexcept :
                vector3D<DirT>{x_dir_comp, y_dir_comp}, pos{x_pos_comp, y_pos_comp} {calc();}
        ray(DirT &&x_pos_comp, DirT &&y_pos_comp, DirT &&x_dir_comp, DirT &&y_dir_comp) noexcept :
                vector3D<DirT>{std::move(x_dir_comp), std::move(y_dir_comp)},
                pos{std::move(x_pos_comp), std::move(y_pos_comp)} {calc();}
        ray(const vector3D<PosT> &origin, DirT &&x_dir_comp, DirT &&y_dir_comp, DirT &&z_dir_comp) noexcept :
                vector3D<DirT>{std::move(x_dir_comp), std::move(y_dir_comp), std::move(z_dir_comp)},
                pos{origin} {calc();}
        void set_pos(const vector3D<PosT> &new_pos) {
            pos = new_pos;
        }
        void set_pos(vector3D<PosT> &&new_pos) {
            pos = std::move(new_pos);
        }
        void set_dir(const vector3D<DirT> &new_dir) { // set_dir() methods need not be called (reassign ray directly)
            *this = new_dir;
            calc();
        }
        void set_dir(vector3D<DirT> &&new_dir) {
            *this = std::move(new_dir);
            calc();
        }
        vector3D<DirT> &get_dir() {
            return *this;
        }
        const vector3D<PosT> &get_pos() const noexcept { // useless method - no need to use
            return pos;
        }
        const LenT &length() const noexcept {
            return l;
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        bool intersects(const body<M, R, T> &bod) {
            /* This method calculates whether the ray intersects with a given body, and if it does, stores the
             * intersection distance in l. The 'b' and 'c' variables below are quadratic formula coefficients - the 'a'
             * coefficient is always equal to 1 (since a = (*this)*(*this)), so has been inlined. */
            this->calc();
            auto b = 2*(pos - bod.curr_pos)*(*this); // operator* between two vectors is treated as scalar product
            auto c = pos*pos - 2*bod.curr_pos*pos + bod.curr_pos*bod.curr_pos - bod.radius*bod.radius;
            auto discriminant = b*b - 4*c; // recall that a = 1
            if (discriminant < 0)
                return false;
            l = minimum((-b + sqrtl(discriminant))/2, (-b - sqrtl(discriminant))/2);
            return true;
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        bool intersects(const body<M, R, T> &&bod) {
            return intersects(bod);
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        bool intersects(const std::vector<std::reference_wrapper<body<M, R, T>>> &bodies) {
            /* This method is similar to the other intersects() overloads, but it will set l to the distance between the
             * ray's origin and the closest body, and return true (if the ray intersects with at least one body). This
             * is intuitive, since a ray in real space will "terminate" once it hits an object, and will (in general)
             * not continue through the object and then intersect with other objects behind it. */
            bool one_intersection = false;
            LenT min_l;
            for (const auto &ref : bodies) {
                if (this->intersects(ref.get())) {
                    if (!one_intersection) {
                        min_l = l;
                        one_intersection = true;
                        continue;
                    }
                    if (l < min_l)
                        min_l = l;
                }
            }
            if (!one_intersection) {
                l = LenT{0};
                return false;
            }
            l = min_l;
            return true;
        }
        ray<PosT, DirT, LenT> &operator=(const vector<DirT> &other) noexcept override {
            vector3D<DirT>::operator=(other);
            this->calc();
            return *this;
        }
        ray<PosT, DirT, LenT> &operator=(const vector3D<DirT> &other) noexcept override {
            vector3D<DirT>::operator=(other);
            this->calc();
            return *this;
        }
        template <typename U> requires (isConvertible<U, DirT>)
        ray<PosT, DirT, LenT> &operator=(const vector3D<U> &other) noexcept {
            vector3D<DirT>::operator=(other);
            this->calc();
            return *this;
        }
        template <typename U> requires (isConvertible<U, DirT>)
        ray<PosT, DirT, LenT> &operator=(const vector3D<U> &&other) noexcept {
            vector3D<DirT>::operator=(other); // no need to use std::move() here
            this->calc();
            return *this;
        }
        ray<PosT, DirT, LenT> &operator=(const vector<DirT> &&other) noexcept override {
            vector3D<DirT>::operator=(std::move(other));
            this->calc();
            return *this;
        }
        ray<PosT, DirT, LenT> &operator=(vector3D<DirT> &&other) noexcept override {
            vector3D<DirT>::operator=(std::move(other));
            this->calc();
            return *this;
        }
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
        friend std::ostream &operator<<(std::ostream &os, const ray<PosU, DirU, LenU> &r);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
        friend std::ostream &operator<<(std::ostream &os, const ray<PosU, DirU, LenU> &&r);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                  isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &&r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator==(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &&r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &&r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &r2);
        template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
                isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
        friend bool operator!=(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &&r2);
    };
    template <isNumWrapper PosT = long double, isNumWrapper DirT = long double, isNumWrapper DistT = long double>
    class camera {
        /* The "camera" class is one which represents a pinhole camera in 3D space with a certain position, orientation
         * (direction it's pointing in) and rotation. A camera is essentially a box with a certain width and height
         * (which are directly proportional to the real image's dimensions) that are "normalised" in order for the
         * shorter of the two to always equal 2. The length of the camera "box" is the image distance. The pinhole lies
         * on one of the sides that are perpendicular to the lengthwise dimension, and the location of the pinhole in
         * simulated 3D space is taken as the location of the camera. On the face opposite the pinhole (aperture) lies
         * the camera's "receptor", which is where the image in simulated 3D space is projected onto. This is done by
         * sending a ray out from every single "pixel" on the receptor through the pinhole, and calculating what it
         * intersects with - although the intersection calculation is not done in the camera class (see astro_scene). */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        static inline const vector3D<DirT> zero{}; // as this is repeatedly used
        static inline const vector3D<DirT> up{DirT{0}, DirT{0}, DirT{1}}; // same with all these
        static inline const vector3D<DirT> down{DirT{0}, DirT{0}, DirT{-1}};
        static inline const vector3D<DirT> right{DirT{1}, DirT{0}, DirT{0}};
        static inline const vector3D<DirT> forwards{DirT{0}, DirT{1}, DirT{0}};
        image_dimensions dims{1000, 1000}; // the real image dimensions in pixels
        vector3D<PosT> pos; // position of the camera
        vector3D<DirT> dir{0, 0, 1}; // direction the camera is pointing in (doesn't have to be a unit vector)
        long double rot = 0; /* this specifies the roll of the camera (see pitch, yaw & roll) in terms of an angle in
        * radians - in the case of the camera's direction being perpendicular to the x-y plane, this angle is with
        * respect to the positive y-axis */
        DistT dist = DistT{2}; /* distance from the camera's "receptor" - i.e., the image projected onto its
        * back - to the camera's aperture (the pinhole) -> greater image distance = greater zoom */
        vector3D<DirT> perp; // vector around which rays are rotated to be in the correct orientation
        long double angle_with_down = 0; // the angle the camera direction makes with the negative z-axis
        long double rotation_correction = 0; // a correction factor needed to due to rays being rotated
        bool points_up = true; // true if the camera is pointing directly up (direction parallel to positive z-axis)
        bool points_down = false;// true if the camera is pointing directly down (direction parallel to negative z-axis)
        long double w_factor = 1; // these factors dictate how the x- and y-components of the ray vectors should be
        long double h_factor = 1; // ... scaled based on the true dimensions of the image
        bool invert_ray = true; /* this parameter dictates whether a ray emanating from the pinhole has its y-component
        * flipped or not - this, essentially, avoids the image having to be flipped later - although for images that
        * have their origin in the top-left corner (unlike BMPs), this should be set to false */
        void calc_all() {
            check_non_zero();
            if (dir.x == DirT{0} && dir.y == DirT{0}) {
                perp.set_x(DirT{1});
                perp.set_y(DirT{0});
                perp.set_z(DirT{0});
                rotation_correction = 0;
                if ((points_down = !(points_up = dir.z > DirT{0})))
                    angle_with_down = 0;
                else
                    angle_with_down = PI;
                goto end;
            }
            perp = cross(down, dir);
            angle_with_down = angle_between(down, dir);
            if (dir.x < 0)
                rotation_correction = angle_between(forwards, dir.xy_projection());
            else if (dir.x > 0)
                rotation_correction = -angle_between(forwards, dir.xy_projection());
            else if (dir.y < 0)
                rotation_correction = PI;
            end:
            if (dims.x >= dims.y) {
                w_factor = ((long double) dims.x)/dims.y;
                h_factor = 1;
            }
            else {
                w_factor = 1;
                h_factor = ((long double) dims.y)/dims.x;
            }
            if (!invert_ray)
                h_factor = -h_factor;
        }
        void check_non_zero() { // the camera cannot, however, not be pointing in any direction
            if (dir == zero)
                throw no_direction("The camera's direction is a zero vector. The camera must point in a direction.\n");
            if (dist == DistT{0})
                throw zero_image_distance("A camera object's image distance cannot be zero, or else simulated bodies"
                                          "directly\nin front of the camera would appear to be infinitely far away.\n");
        }
    public:
        camera() = default; // for a camera at the origin pointing up (taking the z-axis as up)
        explicit camera(const image_dimensions &image_dim) : dims{image_dim} {calc_all();}
        explicit camera(const vector3D<PosT> &position) : pos{position} {calc_all();} // for a camera pointing up
        explicit camera(const vector3D<PosT> &&position) : pos{std::move(position)} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction) :
        pos{position}, dir{direction} {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction) :
        pos{std::move(position)}, dir{std::move(direction)} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction, const long double &rotation,
               const DistT &image_distance) : pos{position}, dir{direction}, rot{rotation}, dist{image_distance}
               {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction, const long double &&rotation,
               const DistT &&image_distance) :
                pos{std::move(position)}, dir{std::move(direction)}, rot{rotation}, dist{std::move(image_distance)}
        {calc_all();}
        camera(const vector3D<PosT> &position, const image_dimensions &image_dim) : pos{position}, dims{image_dim}
        {calc_all();}
        camera(const vector3D<PosT> &&position, const image_dimensions &&image_dim) :
                pos{std::move(position)}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction, const image_dimensions &image_dim) :
                pos{position}, dir{direction}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction, const image_dimensions &&image_dim) :
                pos{std::move(position)}, dir{std::move(direction)}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &position, const vector3D<DirT> &direction, const long double &rotation,
               const DistT &image_distance, const image_dimensions &image_dim) : pos{position}, dir{direction},
               rot{rotation}, dist{image_distance}, dims{image_dim} {calc_all();}
        camera(const vector3D<PosT> &&position, const vector3D<DirT> &&direction, const long double &&rotation,
               const DistT &&image_distance, const image_dimensions &&image_dim) : pos{std::move(position)},
               dir{std::move(direction)}, rot{rotation}, dist{std::move(image_distance)}, dims{image_dim}
        {calc_all();}
        const vector3D<PosT> &position() const noexcept {
            return pos;
        }
        const vector3D<DirT> &direction() const noexcept {
            return dir;
        }
        const long double &rotation() const noexcept {
            return rot;
        }
        const DistT &image_distance() const noexcept {
            return dist;
        }
        const image_dimensions &img_dim() const noexcept {
            return dims;
        }
        camera<PosT, DirT, DistT> &set_position(const vector3D<PosT> &new_position) {
            pos = new_position;
            calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_direction(const vector3D<DirT> &new_direction) {
            dir = new_direction;
            calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_rotation(const long double &new_rotation) {
            rot = new_rotation;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_distance(const DistT &new_image_distance) {
            dist = new_image_distance;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_dimensions(const image_dimensions &new_dimensions) {
            dims = new_dimensions;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_position(vector3D<PosT> &&new_position) {
            pos = std::move(new_position);
            calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_direction(vector3D<DirT> &&new_direction) {
            dir = std::move(new_direction);
            calc_all();
            return *this;
        }
        camera<PosT, DirT, DistT> &set_rotation(const long double &&new_rotation) {
            rot = new_rotation;
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_distance(DistT &&new_image_distance) {
            dist = std::move(new_image_distance);
            return *this;
        }
        camera<PosT, DirT, DistT> &set_image_dimensions(const image_dimensions &&new_dimensions) {
            dims = new_dimensions;
            return *this;
        }
        void recalculate() {
            calc_all();
        }
        camera<PosT, DirT, DistT> &invert_rays(bool invert) {
            invert_ray = invert;
            h_factor = -h_factor;
            return *this;
        }
        long double fovh_rad() const noexcept { // horizontal field of view in radians
            return dims.x <= dims.y ? 2*atanl(1.0l/dist) : 2*atanl(dims.x/(dims.y*dist));
        }
        long double fovv_rad() const noexcept { // vertical field of view in radians
            return dims.y <= dims.x ? 2*atanl(1.0l/dist) : 2*atanl(dims.y/(dims.x*dist));
        }
        long double fovd_rad() const noexcept { // diagonal field of view in radians
            long double x_ratio;
            long double y_ratio;
            if (dims.x > dims.y) {
                x_ratio = ((long double) dims.x)/dims.y;
                y_ratio = 1;
            }
            else {
                x_ratio = 1;
                y_ratio = ((long double) dims.y)/dims.x;
            }
            return 2*atanl(sqrtl(x_ratio*x_ratio + y_ratio*y_ratio)/dist);
        }
        long double fovh_deg() const noexcept { // horizontal field of view in degrees
            return rad_to_deg(fovh_rad());
        }
        long double fovv_deg() const noexcept { // vertical field of view in degrees
            return rad_to_deg(fovv_rad());
        }
        long double fovd_deg() const noexcept { // diagonal field of view in degrees
            return rad_to_deg(fovd_rad());
        }
        template <isNumWrapper LenT = long double>
        ray<PosT, DirT, LenT> get_ray_from_pixel(unsigned int x, unsigned int y) {
            if (x >= dims.x || y >= dims.y)
                throw std::out_of_range("The pixel attempting to be accessed is not within the image dimensions.\n");
            ray<PosT, DirT, LenT> // ray starts out as if it had been shot out of the camera facing directly down
            ray{this->pos,
                ((((long double) x)/dims.x)*2 - 1)*w_factor, ((((long double) y)/dims.y)*2 - 1)*h_factor, -dist};
            ray.rodrigues_rotate(perp, angle_with_down); // rotate to camera direction
            ray.rodrigues_rotate(dir, -rotation_correction - rot); // rotate around camera's z-axis appropriately
            return ray;
        }
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        friend class astro_scene;
    };
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
    class astro_scene : public bmp {
    private:
        image_dimensions dimensions{bmp::width, bmp::height};
        unsigned int num_stars = 4*log(bmp::width*bmp::height);
        long double star_radius = sqrt(bmp::width*bmp::width)/1000.0l; // in pixels
        std::vector<point> star_points;
        std::vector<std::reference_wrapper<body<M, R, T>>> bodies;
        bool rendered = false;
        void create_stars(bool reset) {
            if (!num_stars)
                return;
            sc_clr = colors::white;
            circle c;
            c.set_radius(star_radius);
            std::random_device r;
            if (reset) {
                point p;
                star_points.clear();
                for (unsigned int i = 0; i < num_stars; ++i) {
                    p.x = r() % width;
                    p.y = r() % height;
                    c.set_pos(p);
                    bmp::draw_circle(c);
                    star_points.push_back(p);
                }
                return;
            }
            for (const point &p : star_points) {
                c.set_pos(p);
                bmp::draw_circle(c);
            }
        }
    public:
        astro_scene() :
        bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp")} {}
        explicit astro_scene(const char *source_bmp) :
        bmp{source_bmp,
            std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp")} {}
        astro_scene(unsigned int bmp_width, unsigned int bmp_height) :
        bmp{std::move(get_home_path<String>() + file_sep() + "AstroScene_" + get_date_and_time() + ".bmp"),
            colors::black, bmp_width, bmp_height} {}
        astro_scene(const String &bmp_path, unsigned int bmp_width, unsigned int bmp_height) :
                bmp{bmp_path, colors::black, bmp_width, bmp_height} {}
        astro_scene(String &&bmp_path, color background_color, unsigned int bmp_width, unsigned int bmp_height) :
                bmp{std::move(bmp_path), colors::black, bmp_width, bmp_height} {}
        astro_scene(const astro_scene &other) : bmp{other} {} // not = default, as I will add stuff in body later
        astro_scene(astro_scene &&other) noexcept : bmp{std::move(other)} {}
        bool set_num_stars(unsigned int num) {
            if (num*star_radius*star_radius >= bmp::width*bmp::height)
                return false;
            num_stars = num;
            return true;
        }
        bool set_star_rad(long double rad) {
            if (num_stars*rad*rad >= bmp::width*bmp::height || rad < 0)
                return false;
            star_radius = rad;
            return true;
        }
        bool render(bool reset_star_positions = true) noexcept {
            if (rendered) {
                bmp::sc_clr = colors::black;
                bmp::fill_bg();
            }
            create_stars(reset_star_positions);
            rendered = true;
            return true;
        }
    };
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
    std::ostream &operator<<(std::ostream &os, const ray<PosU, DirU, LenU> &r) {
        return os << '(' << r.pos << ") + " << r.l  << '(' << static_cast<const vector3D<DirU>&>(r) << ')';
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU>
    std::ostream &operator<<(std::ostream &os, const ray<PosU, DirU, LenU> &&r) {
        return os << '(' << r.pos << ") + " << r.l  << '(' << static_cast<const vector3D<DirU>&>(r) << ')';
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator==(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos == r2.pos && r1.l == r2.l &&
               static_cast<const vector3D<DirU>&>(r1) == static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
    template <isNumWrapper PosU, isNumWrapper DirU, isNumWrapper LenU,
              isNumWrapper PosV, isNumWrapper DirV, isNumWrapper LenV>
    bool operator!=(const ray<PosU, DirU, LenU> &&r1, const ray<PosV, DirV, LenV> &&r2) {
        return r1.pos != r2.pos || r1.l != r2.l ||
               static_cast<const vector3D<DirU>&>(r1) != static_cast<const vector3D<DirV>&>(r2);
    }
}
#endif
