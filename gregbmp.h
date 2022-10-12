//
// Created by mario on 10/10/2022.
//

#ifndef GREGBMP_H
#define GREGBMP_H

#include "gregbod.h"

namespace gtd {
#pragma pack(push, 1)
    struct color {
        unsigned char b;
        unsigned char g;
        unsigned char r;
    };
#pragma pack(pop)
    std::ostream &operator<<(std::ostream &os, const color &col);
    namespace colors {
        color black{0, 0, 0};
        color white{255, 255, 255};
        color blue{255, 0, 0};
        color green{0, 255, 0};
        color red{0, 0, 255};
        color pink{180, 105, 255};
        color cerise{99, 49, 222};
        color fuchsia{255, 0, 255};
        color neon_pink{240, 16, 255};
        color pink_orange{128, 152, 248};
        color purple{128, 0, 128};
        color salmon{114, 128, 250};
        color watermelon_pink{131, 115, 227};
        color orange{0, 165, 255};
        color gold{0, 215, 255};
        color yellow{0, 255, 255};
        color lavender{250, 230, 230};
        color indigo{130, 0, 75};
        color violet{238, 130, 238};
        color lime_green{50, 205, 50};
        color forest_green{34, 139, 34};
        color dark_green{0, 100, 0};
        color aqua{255, 255, 0};
        color sky_blue{235, 206, 135};
        color royal_blue{225, 105, 65};
        color navy{128, 0, 0};
        color wheat{179, 222, 245};
        color tan{140, 180, 210};
        color rosy_brown{143, 143, 188};
        color peru{63, 133, 205};
        color chocolate{30, 105, 210};
        color brown{42, 42, 165};
        color maroon{0, 0, 128};
        color snow{250, 250, 255};
        color honey_dew{240, 255, 240};
        color azure{255, 255, 240};
        color ghost_white{255, 248, 248};
        color beige{220, 245, 245};
        color ivory{240, 255, 255};
        color gainsboro{220, 220, 220};
        color silver{192, 192, 192};
        color gray{128, 128, 128};
        color slate_gray{144, 128, 112};
    }
    class zero_size_bmp : public std::invalid_argument {
    public:
        zero_size_bmp() : std::invalid_argument{"The resulting BMP would have a size of zero. Width and height must be "
                                                "non-zero.\n"} {}
    };
    class invalid_bmp_format : public std::exception {
        String msg;
    public:
        invalid_bmp_format() : msg{"BMP is not in the correct format.\n"} {}
        invalid_bmp_format(const char *str) : msg{str} {}
        const char *what() const noexcept override {
            return msg.c_str();
        }
    };
    class bmp {
        static constexpr unsigned int def_width = 1024;
        static constexpr unsigned int def_height = 720;
        unsigned int width = def_width;
        unsigned int height = def_height;
        color **data = nullptr;
        String path;
        color sc_clr = colors::white; // selected colour
        unsigned int stroke_t = 10;
        unsigned char padding[3] = {0, 0, 0};
        unsigned char (*calc_pad)(const unsigned int&) =
                [](const unsigned int &w){return (unsigned char) ((w*3) % 4 == 0 ? 0 : 4 - ((w*3) % 4));};
        void fill_background() {
            data = new color*[height];
            color **ptr = data;
            color *col;
            for (unsigned int j = 0; j < height; ++j, ++ptr) {
                *ptr = new color[width];
                col = *ptr;
                for (unsigned int i = 0; i < width; ++i, ++col) {
                    *col = sc_clr;
                }
            }
        }
        void check_size() {
            if (width == 0 || height == 0) {
                throw zero_size_bmp();
            }
        }
        void dealloc() {
            if (data == nullptr) {
                return;
            }
            color **ptr = data;
            for (unsigned int j = 0; j < height; ++j, ++ptr) {
                delete [] *ptr;
            }
            delete [] data;
            data = nullptr;
        }
        void alloc() {
            dealloc();
            data = new color*[height];
            color **ptr = data;
            for (unsigned int j = 0; j < height; ++j, ++ptr) {
                *ptr = new color[width];
            }
        }
        unsigned char *get_hdr() {
            unsigned int size = file_size();
            static unsigned char header[14] = {'B', 'M', (unsigned char) (size), (unsigned char) (size >> 8),
                                               (unsigned char) (size >> 16), (unsigned char) (size >> 24), 0, 0, 0, 0,
                                               54, 0, 0, 0};
            return header;
        }
        unsigned char *get_info_hdr() {
            static unsigned char info_hdr[40] = {40, 0, 0, 0, (unsigned char) width, (unsigned char) (width >> 8),
                                                 (unsigned char) (width >> 16), (unsigned char) (width >> 24),
                                                 (unsigned char) height, (unsigned char) (height >> 8), (unsigned char)
                                                 (height >> 16), (unsigned char) (height >> 24), 1, 0, 24, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            return info_hdr;
        }
        void read_bmp(const char *source) {
            if (source == nullptr || *source == 0) {
                throw std::invalid_argument("A nullptr or empty string cannot be passed as a source BMP path.\n");
            }
            std::ifstream in(source, std::ios_base::binary | std::ios_base::in);
            if (!in.good()) {
                String msg = "There was an error opening the file \"";
                msg.append_back(source).append_back("\".\n");
                throw std::invalid_argument(msg.c_str());
            }
            unsigned int offset;
            unsigned int info_hdr_size;
            unsigned short bpp;
            in.seekg(10, std::ios_base::beg);
            in.read((char *) &offset, 4);
            in.read((char *) &info_hdr_size, 4);
            in.read((char *) &width, 4);
            in.read((char *) &height, 4);
            in.seekg(2, std::ios_base::cur);
            in.read((char *) &bpp, 2);
            if (bpp != 24) {
                throw invalid_bmp_format("The specified BMP must have a pixel depth of 24 bpp.\n");
            }
            unsigned char pad = calc_pad(width);
            if (!in.good()) {
                bad:
                String msg = "There was an error reading from the file \"";
                msg.append_back(source).append_back("\".\n");
                throw std::invalid_argument(msg.c_str());
            }
            if (offset != 54 || info_hdr_size != 40) {
                throw invalid_bmp_format();
            }
            check_size();
            in.seekg(54, std::ios_base::beg);
            alloc();
            color **ptr = data;
            std::streamsize triple_width = 3*width;
            for (unsigned int j = 0; j < height; ++j) {
                in.read((char *) *ptr++, triple_width);
                in.seekg(pad, std::ios_base::cur);
            }
            if (!in.good())
                goto bad;
            in.close();
        }
        void copy_pixels(const color *const *source) {
            if (source == nullptr) {
                return;
            }
            for (unsigned int j = 0; j < height; ++j) {
                for (unsigned int i = 0; i < width; ++i) {
                    data[j][i] = source[j][i];
                }
            }
        }
    public:
        bmp() : path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_On_" + get_date_and_time() + ".bmp")}
        {fill_background();}
        bmp(const char *source_bmp) :
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_On_" + get_date_and_time() + ".bmp")} {
            read_bmp(source_bmp);
        }
        bmp(unsigned int bmp_width, unsigned int bmp_height) : width{bmp_width}, height{bmp_height},
        sc_clr{colors::black},
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_On_" + get_date_and_time() + ".bmp")} {
            check_size();
            fill_background();
        }
        bmp(color background_color) : sc_clr{background_color},
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_On_" + get_date_and_time() + ".bmp")}
        {fill_background();}
        bmp(color background_color, unsigned int bmp_width, unsigned int bmp_height) : sc_clr{background_color},
        width{bmp_width}, height{bmp_height},
        path{std::move(get_home_path<String>() + file_sep() + "BMP_Image_On_" + get_date_and_time() + ".bmp")} {
            check_size();
            fill_background();
        }
        bmp(const String &bmp_path, color background_color, unsigned int bmp_width, unsigned int bmp_height) :
                path{bmp_path}, sc_clr{background_color}, width{bmp_width}, height{bmp_height} {
            check_size();
            fill_background();
        }
        bmp(String &&bmp_path, color background_color, unsigned int bmp_width, unsigned int bmp_height) :
        path{std::move(bmp_path)}, sc_clr{background_color}, width{bmp_width}, height{bmp_height} {
            check_size();
            fill_background();
        }
        bmp(const bmp &other) : width{other.width}, height{other.height}, path{other.path}, sc_clr{other.sc_clr},
        stroke_t{other.stroke_t} {
            copy_pixels(other.data);
        }
        bmp(bmp &&other) : width{other.width}, height{other.height}, path{std::move(other.path)}, sc_clr{other.sc_clr},
        stroke_t{other.stroke_t} {
            data = other.data;
            other.data = nullptr;
        }
        static bmp create_bmp(const char *source_bmp) {
            return {source_bmp};
        }
        unsigned int get_height() {
            return height;
        }
        unsigned int get_width() {
            return width;
        }
        void set_color(const color &col) {
            sc_clr = col;
        }
        void set_color(const color &&col) {
            sc_clr = col;
        }
        void set_color(unsigned char r, unsigned char g, unsigned char b) {
            sc_clr.r = r;
            sc_clr.g = g;
            sc_clr.b = b;
        }
        bool set_pixel(unsigned int x, unsigned int y, const color &col) {
            if (x >= width || y >= height) {
                return false;
            }
            data[y][x] = col;
            return true;
        }
        bool set_pixel(unsigned int x, unsigned int y, const color &&col) {
            return set_pixel(x, y, col);
        }
        bool set_pixel(unsigned int x, unsigned int y) {
            return set_pixel(x, y, sc_clr);
        }
        color get_color() {
            return sc_clr;
        }
        color get_color(unsigned int x, unsigned int y) {
            if (x >= width || y >= height) {
                return {};
            }
            return data[y][x];
        }
        bool fill_rect(unsigned int x, unsigned int y, unsigned int w, unsigned int h) {
            if (x >= width || y >= height || w == 0 || h == 0) {
                return false;
            }
            if (x + w - 1 >= width) {
                w = width - x;
            }
            if (y + h - 1 >= height) {
                h = height - y;
            }
            unsigned int x_end = x + w;
            unsigned int y_end = y + h;
            for (unsigned int j = y; j < y_end; ++j) {
                for (unsigned int i = x; i < x_end; ++i)
                    data[j][i] = sc_clr;
            }
            return true;
        }
        void fill_bg() { // convenience method
            fill_rect(0, 0, width, height);
        }
        unsigned int get_rgb(unsigned int x, unsigned int y) {
            if (x >= width || y >= height) {
                return -1;
            }
            return (data[y][x].r << 16) + (data[y][x].g << 8) + data[y][x].b;
        }
        void draw_2d_circle(long double x, long double y, long double r) {
            if (r == 0 || x + r < 0 || y + r < 0 || x - r >= width || y - r >= height) {
                return;
            }
            long double diameter = 2*r;
            unsigned int top = roundl(y + r);
            unsigned int bottom = roundl(y - r);
            unsigned int left;
            unsigned int right;
            long double root;
            while (bottom <= top && y < top) {
                root = sqrtl(r*r - (y - bottom)*(y - bottom));
                left = roundl(x - root);
                right = roundl(x + root);
                if (left >= width) {
                    goto end;
                }
                while (left <= right) {
                    data[bottom][left++] = sc_clr;
                }
                end:
                ++bottom;
                if (bottom >= height) {
                    break;
                }
            }
        }
        void make_mono(unsigned char cutoff, color dark_color = colors::black, color bright_color = colors::white) {
            for (unsigned int j = 0; j < height; ++j) {
                for (unsigned int i = 0; i < width; ++i) {
                    if ((data[j][i].r + data[j][i].g + data[j][i].b)/3.0l > cutoff) {
                        data[j][i] = bright_color;
                    }
                    else {
                        data[j][i] = dark_color;
                    }
                }
            }
        }
        void make_grayscale() {
            unsigned char avg;
            for (unsigned int j = 0; j < height; ++j) {
                for (unsigned int i = 0; i < width; ++i) {
                    avg = (data[j][i].r + data[j][i].g + data[j][i].b)/3;
                    data[j][i].r = data[j][i].g = data[j][i].b = avg;
                }
            }
        }
        void clear() {
            dealloc();
        }
        bool set_path(const char *str) {
            if (str == nullptr || *str == 0) {
                return false;
            }
            path = str;
            return true;
        }
        unsigned int file_size() { // has to be always unsigned int, as the BMP format does not allow otherwise
            return 54 + 3*width*height + height*calc_pad(width);
        }
        void reallocate() {
            if (data != nullptr) {
                return;
            }
            alloc();
        }
        bool read(const char *source_bmp) {
            if (source_bmp == nullptr || *source_bmp == 0) {
                return false;
            }
            read_bmp(source_bmp);
        }
        bool write(const char *path_to_bmp) {
            //std::cout << path_to_bmp << std::endl;
            if (path_to_bmp == nullptr || *path_to_bmp == 0) {
                return false;
            }
            std::ofstream out(path_to_bmp, std::ios_base::trunc | std::ios_base::binary);
            if (!out.good()) {
                return false;
            }
            unsigned char pad = calc_pad(width);
            out.write((char *) get_hdr(), 14);
            out.write((char *) get_info_hdr(), 40);
            color **ptr = data;
            std::streamsize triple_width = 3*width;
            for (unsigned int j = 0; j < height; ++j) {
                out.write((char *) *ptr++, triple_width);
                out.write((char *) padding, pad);
            }
            out.close();
            return true;
        }
        bool write() {
            std::cout << path << std::endl;
            return write(path.c_str());
        }
        bmp &operator=(const bmp &other) {
            if (&other == this) {
                return *this;
            }
            bool need_to_alloc = false;
            if (other.data == nullptr || this->width != other.width || this->height != other.height) {
                this->clear();
                need_to_alloc = true;
            }
            this->width = other.width;
            this->height = other.height;
            this->sc_clr = other.sc_clr;
            this->path = other.path;
            this->stroke_t = other.stroke_t;
            if (other.data != nullptr) {
                if (need_to_alloc) {
                    alloc();
                }
                copy_pixels(other.data);
            }
        }
        bmp &operator=(bmp &&other) {
            if (&other == this) { // could happen!
                return *this;
            }
            this->data = other.data;
            other.data = nullptr;
            this->width = other.width;
            this->height = other.height;
            this->sc_clr = other.sc_clr;
            this->path = other.path;
            this->stroke_t = other.stroke_t;
        }
        ~bmp() {
            dealloc();
        }
        friend std::ostream &operator<<(std::ostream &os, const bmp &image);
    };
    std::ostream &operator<<(std::ostream &os, const bmp &image) {
        os << "[gtd::bmp@" << &image << ":width=" << image.width << ",height=" << image.height
           << ",size_allocated=";
        if (image.data == nullptr) {
            return os << "0]";
        }
        return os << image.width*image.height*3 << ']';
    }
    std::ostream &operator<<(std::ostream &os, const color &col) {
        return os << "[gtd::color@" << &col << ":r=" << +col.r << ",g=" << +col.g << ",b=" << +col.b << ']';
    }
}
#endif