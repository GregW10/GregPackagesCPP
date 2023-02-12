#ifndef GREGALG_H
#define GREGALG_H

/* header file for important constants, concept definitions, algorithms (hence the name) and miscellaneous content */

#include <iostream>
#include <vector>
#include <cmath>

#ifdef __PI__
#undef __PI__
#endif

#define __PI__ 3.14159265358979323846264338327950288419716939937510582097494459230

#define MEAN_AVG(a, b) (((a) + (b))/2)

typedef unsigned long long ull_t; // I know there is a typedef in stdint.h (or cstdint), but I like this one!

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
concept isPrintable = requires (const T &a) {
    {std::declval<std::ostream>() << a} -> std::same_as<std::ostream&>;
};
template <typename F, typename T>
concept binaryPredicate = requires (F f, T a, T b) {
    {f(a, b)} -> std::same_as<bool>;
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
    constexpr long double PI = 3.14159265358979323846264338327950288419716939937510582097494459230;
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
        return std::sqrtl(mean_avg(args*args...) - mean*mean);
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
    template<class ForwardIterator, typename T>
    void fill(ForwardIterator begin, ForwardIterator end, const T &value) {
        if (begin == end || begin == end - 1)
            return;
        while (begin != end)
            *begin++ = value;
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