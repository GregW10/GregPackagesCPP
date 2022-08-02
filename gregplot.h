#include <iostream>
#include "gregstring.h"
#include <vector>
#include <cerrno>
#include <cmath>
#include <ctime>

#ifndef errno
extern int errno;
#endif

#pragma pack(push, 1)

#define NUM_PANES 1
#define BPP 24

#define BMP_INFO_SIZE 14
#define BMP_HEADER_SIZE 40

#define DEF_WIDTH 2272
#define DEF_HEIGHT 1704

#define MIN_WIDTH 200
#define MIN_HEIGHT 150

#define DEF_COLOUR 255

#define PADDING(w) ((3*(w)) % 4 == 0 ? 0 : 4 - (3*(w)) % 4)

namespace gtd {
    template <typename T>
    T get_min(const std::vector<T> &vec) {
        if (vec.empty()) {
            return std::numeric_limits<T>::max();
        }
        T min = vec[0];
        for (const T &elem : vec) {
            if (elem < min) {
                min = elem;
            }
        }
        return min;
    }
    template <typename T>
    T get_max(const std::vector<T> &vec) {
        if (vec.empty()) {
            return std::numeric_limits<T>::min();
        }
        T max = vec[0];
        for (const T &elem : vec) {
            if (elem > max) {
                max = elem;
            }
        }
        return max;
    }
    typedef struct {
        unsigned char B;
        unsigned char G;
        unsigned char R;
    } colour;
    template <typename T>
    class plot {
    private:
        typedef struct {
            unsigned short int header;
            unsigned int fileSize;
            unsigned short int reserved;
            unsigned short int reserved2;
            unsigned int offset;
        } BMP_info;
        typedef struct {
            unsigned int header_size;
            unsigned int width;
            unsigned int height;
            unsigned short int num_panes;
            unsigned short int bpp;
            unsigned int compression;
            unsigned int image_size;
            unsigned int width_res;
            unsigned int height_res;
            unsigned int num_colour_palette;
            unsigned int imp_colours;
        } BMP_header;
        BMP_info info;
        BMP_header header;
        gtd::String path;
        colour background;
        u_int width;
        u_int height;
        u_int fileSize;
        colour **image = nullptr;
        std::vector<std::vector<std::pair<T, T>>> plots;
        unsigned char padding[3] = {0, 0, 0};
        colour black = {0, 0, 0};
        colour blue = {255, 0, 0};
        colour green = {0, 255, 0};
        colour red = {0, 0, 255};
        size_t axes_thickness;
        colour axes_colour = {0, 0, 0};
        u_int xtick_thickness;
        u_int xtick_length;
        u_int num_xticks = 10;
        std::vector<T> xtick_positions;
        u_int ytick_thickness;
        u_int ytick_length;
        u_int num_yticks = 8;
        std::vector<T> ytick_positions;
        u_int origin_x;
        u_int origin_y;
        u_int end_x;
        u_int end_y;
        T plotmin_x;
        T plotmin_y;
        T plotmax_x;
        T plotmax_y;
        T min_x;
        T min_y;
        T max_x;
        T max_y;
        bool alloced = false;
        bool xlim_set = false;
        bool ylim_set = false;
        void fill_structs() {
            fileSize = BMP_INFO_SIZE + BMP_HEADER_SIZE + width*height*(BPP/8) + PADDING(width)*height;
            info = {('M' << 8) + 'B', fileSize, 0, 0, BMP_INFO_SIZE + BMP_HEADER_SIZE};
            header = {BMP_HEADER_SIZE, width, height, NUM_PANES, BPP, 0, 0, 0, 0, 0, 0};
        }
        void fill_background() {
            image = (colour **) malloc(sizeof(colour *)*height);
            for (int y = 0; y < height; ++y) {
                *(image + y) = (colour *) malloc(width*(BPP/8));
                for (int x = 0; x < width; ++x) {
                    *(*(image + y) + x) = background;
                }
            }
            alloced = true;
        }
        inline long double gradient(long double x1, long double y1, long double x2, long double y2) {
            return (y2 - y1) / (x2 - x1);
        }
        inline long double intercept(long double m, long double x_point, long double y_point) {
            return y_point - m*x_point;
        }
        inline long double thickness_x_bound(long double x, long double y, long double thickness) {
            return thickness / (2.0*cosl(atanl(y/x)));
        }
        inline long double thickness_y_bound(long double x, long double y, long double thickness) {
            return thickness / (2.0*cosl(atanl(x/y)));
        }
        int draw_line(u_int start_x, u_int start_y, u_int final_x, u_int final_y, u_int thickness,
                      colour clr = {0, 0, 0}) {
            long long x_diff = final_x - start_x;
            long long y_diff = final_y - start_y;
            if (start_x >= width || start_y >= height || (x_diff == 0 && y_diff == 0)) {
                return -1;
            }
            if (final_x >= width) {
                final_x = width - 1;
            }
            if (final_y >= height) {
                final_y = height - 1;
            }
            long double m;
            long double c;
            long double thickness_bound;
            if (abs(x_diff) >= abs(y_diff)) {
                m = gradient(start_x, start_y, final_x, final_y);
                c = intercept(m, start_x, start_y);
                thickness_bound = thickness_x_bound(abs(x_diff), abs(y_diff), thickness);
                u_int thickness_rounded = (u_int) roundl(thickness_bound);
                if (start_x > final_x) {
                    gtd::swap(start_x, final_x);
                }
                u_int y;
                u_int y_w_start;
                u_int y_w_end;
                for (size_t x = start_x; x <= final_x; ++x) {
                    y = (u_int) roundl(m*x + c);
                    y_w_start = thickness_rounded > y ? 0 : y - thickness_rounded;
                    y_w_end = y + thickness_rounded;
                    //while (y_w_start < 0) ++y_w_start;
                    //while (y_w_start >= height) --y_w_start;
                    if (y_w_start >= height) return 1;
                    //while (y_w_end < 0) ++y_w_end;
                    //while (y_w_end >= height) --y_w_end;
                    if (y_w_end >= height) y_w_end = height - 1;
                    for (; y_w_start <= y_w_end; ++y_w_start) image[y_w_start][x] = clr;
                }
                return 0;
            }
            m = gradient(start_y, start_x, final_y, final_x);
            c = intercept(m, start_y, start_x);
            thickness_bound = thickness_y_bound(abs(x_diff), abs(y_diff), thickness);
            u_int thickness_rounded = (u_int) roundl(thickness_bound);
            if (start_y > final_y) {
                gtd::swap(start_y, final_y);
            }
            u_int x;
            u_int x_w_start;
            u_int x_w_end;
            for (size_t y = start_y; y <= final_y; ++y) {
                x = (u_int) roundl(m*y + c);
                // x_w_start = x - thickness_rounded;
                // x_w_end = x + thickness_rounded;
                // while (x_w_start < 0) ++x_w_start;
                // while (x_w_start >= width) --x_w_start;
                // while (x_w_end < 0) ++x_w_end;
                // while (x_w_end >= width) --x_w_end;
                x_w_start = thickness_rounded > x ? 0 : x - thickness_rounded;
                x_w_end = x + thickness_rounded;
                //while (y_w_start < 0) ++y_w_start;
                //while (y_w_start >= height) --y_w_start;
                if (x_w_start >= width) return 1;
                //while (y_w_end < 0) ++y_w_end;
                //while (y_w_end >= height) --y_w_end;
                if (x_w_end >= width) x_w_end = width - 1;
                for (; x_w_start <= x_w_end; ++x_w_start) image[y][x_w_start] = clr;
            }
            return 0;
        }
        //int draw_vertical_line(u_int start_x_coord, u_int start_y_coord, u_int thickness, u_int length) {
        //    if (start_x_coord >= width || start_y_coord >= height) {
        //        return -1;
        //    }
        //    if (start_x_coord + thickness >= width) {
        //        thickness = width - start_x_coord - 1;
        //    }
        //    if (start_y_coord + length >= height);
        //}
        int draw_quad(u_int btm_left_x, u_int btm_left_y, u_int w, u_int h, colour col = {0, 0, 0}) {
            if (btm_left_x < 0 || btm_left_x >= width || btm_left_y < 0 || btm_left_y >= height) {
                return -1;
            }
            if (btm_left_x + w >= width) w = width - btm_left_x - 1;
            if (btm_left_y + h >= height) h = height - btm_left_y - 1;
            for (u_int y = btm_left_y; y <= btm_left_y + h; ++y) {
                for (u_int x = btm_left_x; x <= btm_left_x + w; ++x) image[y][x] = col;
            }
            return 0;
        }
        inline int draw_square(u_int btm_left_x, u_int btm_left_y, u_int side_length, colour col = {0, 0, 0}) {
            int retval = draw_quad(btm_left_x, btm_left_y, side_length, side_length, col);
            return retval;
        }
        inline void create_axes() {
            origin_x = width / 8;
            origin_y = height / 8;
            end_x = width - origin_x;
            end_y = height - origin_y;
            draw_line(origin_x, origin_y, origin_x, end_y, axes_thickness, axes_colour);
            draw_line(origin_x, origin_y, end_x, origin_y, axes_thickness, axes_colour);
            // draw_quad(origin_x, origin_y, axes_thickness, end_y - origin_y, axes_colour);
            // draw_quad(origin_x, origin_y, end_x - origin_x, axes_thickness, axes_colour);
            draw_square(origin_x - roundl(((long double) axes_thickness) / 2), origin_y -
            roundl(((long double) axes_thickness) / 2), axes_thickness / 2, axes_colour);
        }
        void create_x_ticks() {

        }
        size_t get_max_points() {
            size_t retval = 0;
            size_t size;
            for (const std::vector<std::pair<T, T>> &plot : plots) {
                size = plot.size(); // in case of multiple plots, it is faster to only call size() once per iteration
                if (size > retval) {
                    retval = size;
                }
            }
            return retval;
        }
        void set_mins() {
            min_x = plots[0][0].first;
            min_y = plots[0][0].second;
            for (const std::vector<std::pair<T, T>> &plot : plots) {
                for (const std::pair<T, T> &pair : plot) {
                    if (pair.first < min_x) {
                        min_x = pair.first;
                    }
                    if (pair.second < min_y) {
                        min_y = pair.second;
                    }
                }
            }
        }
        void set_maxes() {
            max_x = plots[0][0].first;
            max_y = plots[0][0].second;
            for (const std::vector<std::pair<T, T>> &plot : plots) {
                for (const std::pair<T, T> &pair : plot) {
                    if (pair.first > max_x) {
                        max_x = pair.first;
                    }
                    if (pair.second > max_y) {
                        max_y = pair.second;
                    }
                }
            }
        }
        void write_image(FILE *fp) {
            fwrite(&info, sizeof(BMP_info), 1, fp);
            fwrite(&header, sizeof(BMP_header), 1, fp);
            for (int y = 0; y < height; ++y) {
                fwrite(*(image + y), sizeof(colour), width, fp);
                fwrite(padding, sizeof(unsigned char), PADDING(width), fp);
            }
            for (size_t i = 0; i < height; ++i) {
                free(*(image + i));
            }
            free(image);
            alloced = false;
        }
        inline int check_bmp() {
            if (path.rsubstr(".bmp").empty()) {
                return -1;
            }
            return 0;
        }
        inline void check_T() {
            if (!std::is_integral<T>::value && !std::is_floating_point<T>::value) {
                throw std::invalid_argument("Only primitive data types (i.e., numerical types) can be passed as data "
                                            "for a plot.");
            }
        }
        inline void check_dim() {
            if (width < MIN_WIDTH) {
                width = MIN_WIDTH;
            }
            if (height < MIN_HEIGHT) {
                height = MIN_HEIGHT;
            }
        }
        inline void set_defaults() {
            axes_thickness = width >= height ? height / 200 : width / 200;
            axes_thickness = axes_thickness == 0 ? 1 : axes_thickness;
            ytick_thickness = (xtick_thickness = axes_thickness / 2 == 0 ? 1 : axes_thickness / 2);
            ytick_length = (xtick_length = axes_thickness*2);
        }
    public:
        plot() : path(gtd::get_home_path<gtd::String>()), background{DEF_COLOUR, DEF_COLOUR, DEF_COLOUR},
                 width(DEF_WIDTH),
                 height(DEF_HEIGHT) {
#ifdef _WIN32
            path.append_back("\\Figure1.bmp");
#else
            path.append_back("/Figure1.bmp");
#endif
            check_T();
            set_defaults();
        }

        explicit plot(const char *path_to_fig) :
                path(path_to_fig), background{DEF_COLOUR, DEF_COLOUR, DEF_COLOUR}, width(DEF_WIDTH),
                height(DEF_HEIGHT) {check_T(); set_defaults();}

        plot(const char *path_to_fig, unsigned char R_background, unsigned char G_background,
             unsigned char B_background) :
                path(path_to_fig), background{B_background, G_background, R_background}, width(DEF_WIDTH),
                height(DEF_HEIGHT) {check_T(); set_defaults();}

        plot(const char *path_to_fig, unsigned char R_background, unsigned char G_background,
             unsigned char B_background, u_int fig_width, u_int fig_height) : path(path_to_fig),
             background{B_background, G_background, R_background}, width(fig_width), height(fig_height) {
            check_T();
            check_dim();
            set_defaults();
        }

        plot(const char *path_to_fig, colour background_colour) : path(path_to_fig), background(background_colour),
        width(DEF_WIDTH), height(DEF_HEIGHT) {check_T(); set_defaults();}

        plot(const char *path_to_fig, colour background_colour, u_int fig_width, u_int fig_height) : path(path_to_fig),
        background(background_colour), width(fig_width), height(fig_height) {check_T(); check_dim(); set_defaults();}

        int add_plot(std::vector<T> x, std::vector<T> y) {
            if (x.size() == 0 || y.size() == 0) {
                return -1;
            }
            size_t size = x.size() < y.size() ? x.size() : y.size();
            std::vector<std::pair<T, T>> plot(size);
            size_t index = 0;
            for (std::pair<T, T> &p : plot) {
                p = {x[index], y[index]};
                ++index;
            }
            plots.push_back(plot);
            return 0;
        }

        void display() const noexcept {
            size_t num = 1;
            for (const auto &plot : plots) {
                std::cout << "Plot " << num++ << ":\n";
                for (const auto &pair : plot) {
                    std::cout << "x: " << pair.first << ", y: " << pair.second << '\n';
                }
                std::cout << "\n";
            }
            std::cout << std::endl;
        }

        void clear_plots() noexcept {
            plots.clear();
        }

        int set_num_xticks(u_int number_of_xticks) {
            if (number_of_xticks*xtick_thickness >= end_x - origin_x) {
                return -1;
            }
            num_xticks = number_of_xticks;
            return 0;
        }

        int set_num_yticks(u_int number_of_yticks) {
            if (number_of_yticks*ytick_thickness >= end_y - origin_y) {
                return -1;
            }
            num_yticks = number_of_yticks;
            return 0;
        }

        int set_xtick_thickness(u_int thickness) {
            if (thickness > axes_thickness) {
                return -1;
            }
            xtick_thickness = thickness;
            return 0;
        }

        int set_ytick_thickness(u_int thickness) {
            if (thickness > axes_thickness) {
                return -1;
            }
            ytick_thickness = thickness;
            return 0;
        }

        int set_xtick_positions(const std::vector<T> &positions) {
            size_t size = positions.size();
            if (size*xtick_thickness >= end_x - origin_x) {
                return -1;
            }
            num_xticks = size;
            xtick_positions = positions;
            return 0;
        }

        int set_xtick_positions(const T *positions) { // must be terminated by the largest possible value of T
            T max = std::numeric_limits<T>::max();
            while (*positions < max) {
                xtick_positions.push_back(*positions++);
            }
            size_t size = xtick_positions.size();
            if (size*xtick_thickness >= end_x - origin_x) {
                return -1;
            }
            num_xticks = size;
            return 0;
        }

        int set_ytick_positions(const std::vector<T> &positions) {
            size_t size = positions.size();
            if (size*ytick_thickness >= end_y - origin_y) {
                return -1;
            }
            num_yticks = size;
            ytick_positions = positions;
            return 0;
        }

        int set_ytick_positions(const T *positions) { // must be terminated by the largest possible value of T
            T max = std::numeric_limits<T>::max();
            while (*positions < max) {
                ytick_positions.push_back(*positions++);
            }
            size_t size = ytick_positions.size();
            if (size*ytick_thickness >= end_y - origin_y) {
                return -1;
            }
            num_yticks = size;
            return 0;
        }

        int set_xlim(T lower_bound, T upper_bound) { // in terms of T units
            if (lower_bound >= upper_bound) {
                return -1;
            }
            plotmin_x = lower_bound;
            plotmax_x = upper_bound;
            xlim_set = true;
            return 0;
        }

        int set_ylim(T lower_bound, T upper_bound) { // in terms of T units
            if (lower_bound >= upper_bound) {
                return -1;
            }
            plotmin_y = lower_bound;
            plotmax_y = upper_bound;
            ylim_set = true;
            return 0;
        }

        int set_axes_thickness(u_int thickness) { // in terms of pixels
            if (thickness > height/20 || thickness > width/20) {
                return -1;
            }
            axes_thickness = thickness;
            return 0;
        }

        void set_axes_colour(colour col) {
            axes_colour = col;
        }

        int gen_plot() {
            if (check_bmp() == -1 || plots.empty()) {
                errno = EINVAL;
                return -1;
            }
            FILE *fp = fopen(path.c_str(), "wb+");
            if (fp == nullptr) {
                errno = EIO;
                return -1;
            }
            check_dim();
            fill_structs();
            fill_background();
            create_axes();
            create_x_ticks();
            draw_line(500, 500, 5500, 5500, 5, blue);
            //draw_square(100, 100.0F, 500, {50, 50, 255});
            // int retval = draw_quad(1000, 1000, 6000, 6200, green);
            write_image(fp);
            fclose(fp);
            return 0;
        }

        void *operator[](const char *element) {
            if (element == nullptr || *element == 0) {
                return nullptr;
            }
            char elem[11];
            strcpy_c(elem, element);
            string_upper(elem);
            if (strcmp_c(elem, "PATH") == 0) {
                return &(path[0]);
            } else if (strcmp_c(elem, "BACKGROUND") == 0) {
                return &background;
            } else if (strcmp_c(elem, "WIDTH") == 0) {
                return &width;
            } else if (strcmp_c(elem, "HEIGHT") == 0) {
                return &height;
            }
            return nullptr;
        }
        ~plot() {
            if (alloced) {
                for (size_t i = 0; i < height; ++i) {
                    free(*(image + i));
                }
                free(image);
                alloced = false;
            }
        }
    };
    std::ostream &operator<<(std::ostream &out, const colour &col) {
        if (!out.good()) {
            return out;
        }
        out << "Blue: " << (short int) col.B << ", Green: " << (short int) col.G << ", Red: " << (short int) col.R;
        return out;
    }
    namespace colours {
        colour black{0, 0, 0};
        colour white{255, 255, 255};
        colour blue{255, 0, 0};
        colour green{0, 255, 0};
        colour red{0, 0, 255};
        colour pink{180, 105, 255};
        colour cerise{99, 49, 222};
        colour fuchsia{255, 0, 255};
        colour neon_pink{240, 16, 255};
        colour pink_orange{128, 152, 248};
        colour purple{128, 0, 128};
        colour salmon{114, 128, 250};
        colour watermelon_pink{131, 115, 227};
        colour orange{0, 165, 255};
        colour gold{0, 215, 255};
        colour yellow{0, 255, 255};
        colour lavender{250, 230, 230};
        colour indigo{130, 0, 75};
        colour violet{238, 130, 238};
        colour lime_green{50, 205, 50};
        colour forest_green{34, 139, 34};
        colour dark_green{0, 100, 0};
        colour aqua{255, 255, 0};
        colour sky_blue{235, 206, 135};
        colour royal_blue{225, 105, 65};
        colour navy{128, 0, 0};
        colour wheat{179, 222, 245};
        colour tan{140, 180, 210};
        colour rosy_brown{143, 143, 188};
        colour peru{63, 133, 205};
        colour chocolate{30, 105, 210};
        colour brown{42, 42, 165};
        colour maroon{0, 0, 128};
        colour snow{250, 250, 255};
        colour honey_dew{240, 255, 240};
        colour azure{255, 255, 240};
        colour ghost_white{255, 248, 248};
        colour beige{220, 245, 245};
        colour ivory{240, 255, 255};
        colour gainsboro{220, 220, 220};
        colour silver{192, 192, 192};
        colour gray{128, 128, 128};
        colour slate_gray{144, 128, 112};
    }
}
#pragma pack(pop)