#include <iostream>
#include "gregstring.h"
#include <vector>
#include <cerrno>
#include <cmath>
#include <ctime>

#ifndef errno
extern int errno;
#endif

#pragma pack(1)

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
        u_int tick_thickness;
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
        }
        long double gradient(long double x1, long double y1, long double x2, long double y2) {
            return (y2 - y1) / (x2 - x1);
        }
        long double intercept(long double m, long double x_point, long double y_point) {
            return y_point - m*x_point;
        }
        long double thickness_x_bound(long double x, long double y, long double thickness) {
            long double theta = atan(y/x);
            return thickness / (2*cos(theta));
        }
        long double thickness_y_bound(long double x, long double y, long double thickness) {
            long double theta = atan(x/y);
            return thickness / (2*cos(theta));
        }
        int draw_line(u_int start_x, u_int start_y, u_int end_x, u_int end_y, u_int thickness,
                      colour clr = {0, 0, 0}) {
            long int x_diff = end_x - start_x;
            long int y_diff = end_y - start_y;
            if (start_x > width - 1 || end_x > width - 1 || start_y > height - 1 || end_y > height - 1 ||
                (x_diff == 0 && y_diff == 0)) {
                return -1;
            }
            long double m;
            long double c;
            long double thickness_bound;
            if (abs(x_diff) >= abs(y_diff)) {
                m = gradient(start_x, start_y, end_x, end_y);
                c = intercept(m, start_x, start_y);
                thickness_bound = thickness_x_bound(abs(x_diff), abs(y_diff), thickness);
                u_int thickness_rounded = (u_int) roundl(thickness_bound);
                if (start_x > end_x) {
                    gtd::swap(start_x, end_x);
                }
                u_int y;
                u_int y_w_start;
                u_int y_w_end;
                for (size_t x = start_x; x <= end_x; ++x) {
                    y = (u_int) roundl(m*x + c);
                    y_w_start = y - thickness_rounded;
                    y_w_end = y + thickness_rounded;
                    while (y_w_start < 0) {
                        ++y_w_start;
                    }
                    while (y_w_start > height - 1) {
                        --y_w_start;
                    }
                    while (y_w_end < 0) {
                        ++y_w_end;
                    }
                    while (y_w_end > height - 1) {
                        --y_w_end;
                    }
                    for (; y_w_start <= y_w_end; ++y_w_start) {
                        image[y_w_start][x] = clr;
                    }
                }
                return 0;
            }
            m = gradient(start_y, start_x, end_y, end_x);
            c = intercept(m, start_y, start_x);
            thickness_bound = thickness_y_bound(abs(x_diff), abs(y_diff), thickness);
            u_int thickness_rounded = (u_int) roundl(thickness_bound);
            if (start_y > end_y) {
                gtd::swap(start_y, end_y);
            }
            u_int x;
            u_int x_w_start;
            u_int x_w_end;
            for (size_t y = start_y; y <= end_y; ++y) {
                x = (u_int) roundl(m*y + c);
                x_w_start = x - thickness_rounded;
                x_w_end = x + thickness_rounded;
                while (x_w_start < 0) {
                    ++x_w_start;
                }
                while (x_w_start > width - 1) {
                    --x_w_start;
                }
                while (x_w_end < 0) {
                    ++x_w_end;
                }
                while (x_w_end > width - 1) {
                    --x_w_end;
                }
                for (; x_w_start <= x_w_end; ++x_w_start) {
                    image[y][x_w_start] = clr;
                }
            }
            return 0;
        }
        int draw_quad(u_int btm_left_x, u_int btm_left_y, u_int w, u_int h, colour col = {0, 0, 0}) {
            if (btm_left_x < 0 || btm_left_x > width - 1 || btm_left_y < 0 || btm_left_y > height - 1) {
                return -1;
            }
            if (btm_left_x + w > width - 1) {
                w = width - btm_left_x - 1;
            }
            if (btm_left_y + h > height - 1) {
                h = height - btm_left_y - 1;
            }
            for (u_int y = btm_left_y; y <= btm_left_y + h; ++y) {
                for (u_int x = btm_left_x; x <= btm_left_x + w; ++x) {
                    image[y][x] = col;
                }
            }
            return 0;
        }
        int draw_square(u_int btm_left_x, u_int btm_left_y, u_int side_length, colour col = {0, 0, 0}) {
            int retval = draw_quad(btm_left_x, btm_left_y, side_length, side_length, col);
            return retval;
        }
        inline void create_axes() {
            u_int origin_x = width / 8;
            u_int origin_y = height / 8;
            u_int end_x = width - origin_x;
            u_int end_y = height - origin_y;
            draw_line(origin_x, origin_y, origin_x, end_y, axes_thickness);
            draw_line(origin_x, origin_y, end_x, origin_y, axes_thickness);
            draw_square(origin_x - roundl(((long double) axes_thickness) / 2), origin_y -
            roundl(((long double) axes_thickness) / 2), axes_thickness / 2);
        }
        void create_x_ticks() {
            u_int no_x_ticks;

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
        T *get_mins() {
            T x_min = plots[0][0].first;
            T y_min = plots[0][0].second;
            for (const std::vector<std::pair<T, T>> &plot : plots) {
                for (const std::pair<T, T> &pair : plot) {
                    if (pair.first < x_min) {
                        x_min = pair.first;
                    }
                    if (pair.second < y_min) {
                        y_min = pair.second;
                    }
                }
            }
            T *mins = (T *) malloc(2*sizeof(T));
            *mins = x_min;
            *(mins + 1) = y_min;
            return mins;
        }
        T *get_maxes() {
            T x_max = plots[0][0].first;
            T y_max = plots[0][0].second;
            for (const std::vector<std::pair<T, T>> &plot : plots) {
                for (const std::pair<T, T> &pair : plot) {
                    if (pair.first > x_max) {
                        x_max = pair.first;
                    }
                    if (pair.second > y_max) {
                        y_max = pair.second;
                    }
                }
            }
            T *maxes = (T *) malloc(2*sizeof(T));
            *maxes = x_max;
            *(maxes + 1) = y_max;
            return maxes;
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
        }

        explicit plot(const char *path_to_fig) :
                path(path_to_fig), background{DEF_COLOUR, DEF_COLOUR, DEF_COLOUR}, width(DEF_WIDTH),
                height(DEF_HEIGHT) {check_T();}

        plot(const char *path_to_fig, unsigned char R_background, unsigned char G_background,
             unsigned char B_background) :
                path(path_to_fig), background{B_background, G_background, R_background}, width(DEF_WIDTH),
                height(DEF_HEIGHT) {check_T();}

        plot(const char *path_to_fig, unsigned char R_background, unsigned char G_background,
             unsigned char B_background, u_int fig_width, u_int fig_height) : path(path_to_fig),
             background{B_background, G_background, R_background}, width(fig_width), height(fig_height) {
            check_T();
            check_dim();
        }

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
            for (const auto &plot : plots) {
                for (const auto &pair : plot) {
                    std::cout << pair.first << " " << pair.second << '\n';
                }
                std::cout << std::endl;
            }
        }

        void clear_plots() noexcept {
            plots.clear();
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
            axes_thickness = width >= height ? height / 200 : width / 200;
            axes_thickness = axes_thickness == 0 ? 1 : axes_thickness;
            tick_thickness = axes_thickness / 2 == 0 ? 1 : axes_thickness / 2;
            fill_structs();
            fill_background();
            create_axes();
            create_x_ticks();
            draw_line(500, 500, 550, 550, 5, blue);
            draw_square(100, 100.0F, 500, {50, 50, 255});
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
        ~plot() = default;
    };
    std::ostream &operator<<(std::ostream &out, const colour &col) {
        if (!out.good()) {
            return out;
        }
        out << "Blue: " << (short int) col.B << ", Green: " << (short int) col.G << ", Red: " << (short int) col.R;
        return out;
    }
}