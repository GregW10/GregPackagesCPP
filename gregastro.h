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
    template <isNumWrapper PosT, isNumWrapper DirT, isNumWrapper ApT> // ApT for aperture type
    class camera {
        vector3D<PosT> pos; // position of the camera
        vector3D<DirT> dir{1, 0, 0}; // direction the camera is pointing in (doesn't have to be a unit vector)
        long double rotation = 0; /* this specifies the roll of the camera (see pitch, yaw & roll) in terms of an angle
        * in radians - in the case of the camera's direction being perpendicular to the x-y plane, this angle is with
        * respect to the positive y-axis */
        ApT ap_w{0}; // aperture width
        ApT ap_h{0}; // aperture height
        long double zoom = 0.5; // between 0 and 1 - measure of opening angle
        void check_non_zero_dir() { // the camera cannot, however, not be pointing in any direction
            DirT zero{0};
            if (dir[0] == zero && dir[1] == zero && dir[2] == zero)
                throw no_direction("The camera's direction is a zero vector. The camera must point in a direction.\n");
        }
    public:
        camera() = default; // for a camera at the origin pointing up (taking the z-axis as up)
        explicit camera(const vector3D<long double> &position) : pos{position} {} // for a camera pointing up
        camera(const vector3D<long double> &position, const vector3D<long double> &direction) :
        pos{position}, dir{direction} {check_non_zero_dir();}
        camera(const vector3D<long double> &&position, const vector3D<long double> &&direction) :
        pos{position}, dir{direction} {check_non_zero_dir();}
        template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
        friend class astro_scene;
    };
    template <isNumWrapper M, isNumWrapper R, isNumWrapper T>
    class astro_scene : public bmp {
    private:
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
}
#endif
