#ifndef GREGALG_H
#define GREGALG_H

/* header file for important constants, concept definitions, algorithms (hence the name) and miscellaneous content */

#ifndef __cplusplus
#error "The gregalg.hpp header file is a C++ header file only."
#endif

#include <cstdint>

#ifndef UINT64_MAX
#error "64-bit (fixed width) integral data type not available. Compilation failed."
#endif

#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <cinttypes>
#include <climits>

#ifdef __PI__
#undef __PI__
#endif

#define __PI__ 3.14159265358979323846264338327950288419716939937510582097494459230l

#define BILLION 1'000'000'000.0l

#define MEAN_AVG(a, b) (((a) + (b))/2)

typedef unsigned long long ull_t;

template <typename From, typename To>
concept isConvertible = std::is_convertible<From, To>::value;
template <typename T>
concept isIntegral = std::is_integral<T>::value;
// isNumWrapper is also used for the vector, vector2D and vector3D classes
template <typename T> // any primitive type, or any class acting as a numerical wrapper type which implements...
concept isNumWrapper = requires (T val, T val2, size_t l, long double f, std::ostream out) {
    T{1}; // ...a constructor which accepts numerical types
    T{T{}}; // a constructor which accepts another type T, and a default (empty) constructor (should ideally construct
    // an object corresponding to the number zero, or else numerous methods within matrix<T> do not make sense)
    {val + val2} -> isConvertible<T>; // must be overloaded for addition onto itself
    {val + l} -> isConvertible<T>; // must be overloaded for addition onto integer type
    {l + val} -> isConvertible<T>; // ... and vice versa
    {val += val2} -> isConvertible<T>;
    {val += l} -> isConvertible<T>;
    {val - val2} -> isConvertible<T>; // same for subtraction
    {val - l} -> isConvertible<T>;
    {l - val} -> isConvertible<T>;
    {val -= val2} -> isConvertible<T>;
    {val -= l} -> isConvertible<T>;
    {val*val2} -> isConvertible<T>; // multiplication
    {val*l} -> isConvertible<T>;
    {l*val} -> isConvertible<T>;
    {val *= val2} -> isConvertible<T>;
    {val *= l} -> isConvertible<T>;
    {val/val2} -> isConvertible<T>; // same for division
    {val/l} -> isConvertible<T>;
    {l/val} -> isConvertible<T>;
    {val /= val2} -> isConvertible<T>;
    {val /= l} -> isConvertible<T>;
    {val == l} -> std::same_as<bool>; // must have the comparison operators overloaded
    {val != l} -> std::same_as<bool>;
    {val + f} -> isConvertible<T>; // must be overloaded for addition onto float type
    {f + val} -> isConvertible<T>; // ... and vice versa
    {val += f} -> isConvertible<T>;
    {val - f} -> isConvertible<T>;
    {f - val} -> isConvertible<T>;
    {val -= f} -> isConvertible<T>;
    {val*f} -> isConvertible<T>;
    {f*val} -> isConvertible<T>;
    {val *= f} -> isConvertible<T>;
    {val/f} -> isConvertible<T>;
    {f/val} -> isConvertible<T>;
    {val /= f} -> isConvertible<T>;
    {val == f} -> std::same_as<bool>; // must have the comparison operators overloaded
    {val != f} -> std::same_as<bool>;
    {val == val2} -> std::same_as<bool>;
    {val != val2} -> std::same_as<bool>;
    {out << val}; // must have the insertion operator overloaded for outputting to a std::ostream object
}; // see isIntegralNumWrapper concept in gregvec.h for modulo requirements
template <typename T>
concept isPrintable = requires (std::ostream &os, const T &a) {
    {os << a} -> std::same_as<std::ostream&>;
    // {std::move(os) << a} -> std::same_as<std::ostream&>;
};
template <typename F, typename T>
concept binaryPredicate = requires (F f, T a, T b) {
    {f(a, b)} -> std::same_as<bool>;
};
template <typename F, typename T>
concept unaryPredicate = requires (F f, T a) {
    {f(a)} -> std::same_as<bool>;
};
template <typename IT>
concept forwardIterator = requires (IT it, IT other) {
    {it++} -> std::same_as<IT>;
    {++it} -> std::same_as<IT&>;
    {*it};
    {*it = *other};
    {it != other} -> std::same_as<bool>;
};
namespace gtd {
    constexpr long double PI = 3.14159265358979323846264338327950288419716939937510582097494459230l;
    template<typename T>
    inline void swap(T &A, T &B) {
        T C{A};
        A = B;
        B = C;
    }
    template<typename T>
    inline void swap(T *A, T *B) {
        T C{*A};
        *A = *B;
        *B = C;
    }
    template <typename T> requires requires (T a, T b) {{a <= b} -> std::same_as<bool>;}
    inline bool le(const T &a, const T &b) {
        return a <= b;
    }
    // the below 3 functions emulate 'zip()' from Python (I really like what I've done here if I may say so myself)
    template <typename... pack> // first: creates copies of all the elements in the vectors passed
    std::vector<std::tuple<pack...>> zip_cpy(const std::vector<pack>&... vectors) {
        if constexpr (!sizeof...(pack)) // to reject branches at compile time
            return std::vector<std::tuple<>>();
        else {
            std::vector<std::tuple<pack...>> zipped;
            std::vector<size_t> sizes = {(vectors.size())...};
            size_t size = sizes[0];
            if (!std::all_of(sizes.begin(), sizes.end(), [&size](size_t s) { return s == size; }))
                return zipped;
            for (size_t i = 0; i < size; ++i)
                zipped.push_back(std::tuple<pack...>(vectors[i]...));
            return zipped;
        }
    }
    template <typename... pack> // second: with references
    std::vector<std::tuple<pack&...>> zip_ref(std::vector<pack>&... vectors) {
        if constexpr (!sizeof...(pack))
            return std::vector<std::tuple<>>();
        else {
            std::vector<std::tuple<pack&...>> zipped;
            std::vector<size_t> sizes = {(vectors.size())...};
            size_t size = sizes[0];
            if (!std::all_of(sizes.begin(), sizes.end(), [&size](size_t s) { return s == size; }))
                return zipped;
            for (size_t i = 0; i < size; ++i)
                zipped.push_back(std::move(std::tie(vectors[i]...)));
            return zipped;
        }
    }
    template <typename... pack> // third: with const references
    std::vector<std::tuple<const pack&...>> zip_cref(const std::vector<pack>&... vectors) {
        if constexpr (!sizeof...(pack))
            return std::vector<std::tuple<>>();
        else {
            std::vector<std::tuple<const pack&...>> zipped;
            std::vector<size_t> sizes = {(vectors.size())...};
            size_t size = sizes[0];
            if (!std::all_of(sizes.begin(), sizes.end(), [&size](size_t s) {return s == size;}))
                return zipped;
            for (size_t i = 0; i < size; ++i)
                zipped.push_back(std::tie(std::as_const(vectors[i])...));
            return zipped;
        }
    }
    template <forwardIterator IT>
    auto mean_avg(IT begin, IT end) -> typename std::remove_reference<decltype(*(std::declval<IT>()))>::type {
        if (begin == end)
            return typename std::remove_reference<decltype(*(std::declval<IT>()))>::type{};
        typename std::remove_reference<decltype(*(std::declval<IT>()))>::type total{0};
        ull_t count = 0;
        while (begin != end) {
            total += *begin++;
            ++count;
        }
        return total/count;
    }
    template <forwardIterator IT, isNumWrapper OUT>
    void mean_avg(IT begin, IT end, OUT &out) {
        if (begin == end)
            return;
        out = 0;
        ull_t count = 0;
        while (begin != end) {
            out += *begin++;
            ++count;
        }
        out /= count;
    }
    /* lots of compile-time magic below */
    template <isNumWrapper ...Args>
    constexpr inline auto mean_avg(const Args& ...args) {
        static_assert(sizeof...(args), "mean_avg cannot be called with zero arguments.\n");
        return (args + ...)/sizeof...(args);
    }
    template <isNumWrapper ...Args>
    constexpr inline auto sd(const Args& ...args) {
        static_assert(sizeof...(args), "sd cannot be called with zero arguments.\n");
        auto &&mean = mean_avg(args...);
        /* Although I would prefer using the mathematical functions within the 'std' namespace (such as std::sqrtl),
         * there exist various standard non-conforming GCC (g++) versions (and other compilers) that do not have certain
         * versions of the mathematical functions in the 'std' namespace, and only in the global namespace. */
        return sqrtl(mean_avg(args*args...) - mean*mean);
    }
    template <typename T, binaryPredicate<const T&> F>
    requires requires (T first, T second) {{first < second} -> std::same_as<bool>;}
    void bubble_sort(std::vector<T> &array, const F &compare = [](const T& t1, const T& t2){return t1 <= t2;}) {
        using vec_size_t = typename std::vector<T>::size_type;
        vec_size_t size = array.size();
        if (!size || size == 1)
            return;
        T *data = array.data();
        T *first;
        T *second;
        T *end = data + size;
        vec_size_t count;
        bool sorted = false;
        while (!sorted) {
            count = 1;
            first = data;
            second = first + 1;
            --end; // discard last element from consideration after each sort
            for (; second <= end; ++first, ++second) {
                if (!compare(*first, *second)) {
                    swap(first, second);
                    continue;
                }
                ++count;
            }
            sorted = (count == size--);
        }
    }
    template <typename T, binaryPredicate<const T&> F = bool (*)(const T&, const T&)>
    requires requires (T first, T second) {{first < second} -> std::same_as<bool>;}
    void n_sort(std::vector<T> &array, const F &compare = le) {
        using vec_size_t = typename std::vector<T>::size_type;
        vec_size_t size = array.size();
        if (!size || size == 1)
            return;
        std::vector<T> sorted(size);
        T *a_ptr;
        T *s_ptr = sorted.data();
        T extreme;
        vec_size_t extreme_index;
        vec_size_t outer = 0;
        vec_size_t inner;
        vec_size_t new_size = size;
        while (outer++ < size) {
            a_ptr = array.data();
            extreme = *a_ptr++;
            extreme_index = 0;
            for (inner = 1; inner < new_size; ++inner, ++a_ptr) {
                if (!compare(extreme, *a_ptr)) {
                    extreme = *a_ptr;
                    extreme_index = inner;
                }
            }
            *s_ptr++ = extreme;
            array.erase(array.begin() + extreme_index);
            --new_size;
        }
        array.swap(sorted);
    }
    template <typename T>
    void quick_sort(std::vector<T> &array) {
        using vec_size_t = typename std::vector<T>::size_type;
        vec_size_t size = array.size();
        if (!size || size == 1)
            return;
        // remains to be seen - might have to write a pivot or sorter class
    }
    template<class BiDirectionalIterator>
    void reverse(BiDirectionalIterator begin, BiDirectionalIterator end) {
        end--;
        if (begin == end)
            return;
        while (begin < end)
            swap(*begin++, *end--);
    }
    template<forwardIterator fIT, typename T>
    void fill(fIT begin, fIT end, const T &value) {
        if (begin == end || begin == end - 1)
            return;
        while (begin != end)
            *begin++ = value;
    }
    /* the below function, simple as it is, I find highly useful as it performs a byte-wise swap, and can, therefore,
     * be used on objects without touching any of their copy constructors */
    void swap(void *one, void *two, size_t size) {
        if (!one || !two || !size)
            return;
        char temp;
        char *first = (char *) one;
        char *second = (char *) two;
        while (size --> 0) {
            temp = *first;
            *first++ = *second;
            *second++ = temp;
        }
    }
    void move(void *dest, void *source, size_t size) {
        if (!dest || !source || !size)
            return;
        char *dst = (char *) dest;
        char *src = (char *) source;
        while (size --> 0) {
            *dst++ = *src;
            *src++ = 0;
        }
    }
    void copy(void *dest, const void *source, size_t size) {
        if (!dest || !source || !size)
            return;
        char *dst = (char *) dest;
        char *src = (char *) source;
        while (size --> 0) *dst++ = *src++;
    }
    inline long double rad_to_deg(const long double &radians) {
        return (radians*180)/PI;
    }
    inline long double deg_to_rad(const long double &degrees) {
        return (degrees*PI)/180;
    }
    inline long double rad_to_deg(const long double &&radians) {
        return (radians*180)/PI;
    }
    inline long double deg_to_rad(const long double &&degrees) {
        return (degrees*PI)/180;
    }
    template <isNumWrapper T>
    auto sqrt(const T& val, const T& epsilon = 0.00000001l) { // brute force
        if (val < 0)
            throw std::invalid_argument{"Error: cannot compute square root of a negative number.\n"};
        // if (val == 0)
        //     return 0.0l;
        T guess = val < 10 ? val/2 : val/10;
        T delta = guess/8;
        T prev_guess;
        T val_guess;
        T diff;
        bool is_under = false;
        do {
            prev_guess = guess;
            if ((val_guess = guess*guess) < val) {
                guess += delta;
                if (!is_under)
                    delta /= 8;
                is_under = true;
            }
            else {
                guess -= delta;
                if (is_under)
                    delta /= 8;
                is_under = false;
            }
            diff = val_guess < val ? val - val_guess : val_guess - val;
        } while (diff > epsilon);
        return guess;
    }
}
template <isPrintable T> // must be defined in the global namespace due to argument-dependent-lookup
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    os << '[';
    for (const T &val : vec)
        os << val << ", ";
    return os << "\b\b]";
}
#endif
