//
// Created by mario on 01/09/2022.
//

#ifndef GREGVEC_H
#define GREGVEC_H

#include "gregstring.h"
#include "gregmat.h"

#include <iostream>
#include <type_traits>
#include <utility>

namespace gtd {
    template <isFund T>
    class vector {
    protected:
        static String repr; // only one gtd::String for all vectors to minimise dynamic memory allocation
        typedef T vector_type;
    public:
        vector() = default;
        virtual const char *c_str(unsigned char f_p_dec_places = 5) const = 0;
        virtual long double magnitude() const noexcept = 0;
        virtual const T &operator[](unsigned char index) const noexcept = 0;
        virtual vector<T> &operator++() noexcept = 0; //can only declare reference-returning func. for an abstract class
        virtual vector<T> &operator--() noexcept = 0;
        virtual vector<T> &operator=(const vector<T> &other) noexcept = 0; // copies other
        template <isFund U>
        friend class vector;
    };
    template <isFund T>
    class vector3D;
    template <isFund T>
    class vector2D : public vector<T> {
    protected:
        T x = 0;
        T y = 0;
    public:
        vector2D() = default;
        template <isConvertible<T> U>
        explicit vector2D(const vector2D<U> &other) noexcept : x{other.x}, y{other.y} {}
        // vector2D(T &&x_component, T &&y_component) noexcept : x{x_component}, y{y_component} {}
        //vector2D(T &x_component, T &y_component) noexcept : x{x_component}, y{y_component} {}
        vector2D(T x_component, T y_component) noexcept : x{x_component}, y{y_component} {}
        const char *c_str(unsigned char f_p_dec_places = 5) const override {
            vector<T>::repr.erase_chars();
            vector<T>::repr.append_back(x, f_p_dec_places).append_back("i + ").append_back(y, f_p_dec_places).
            push_back('j');
            return vector<T>::repr.c_str();
        }
        virtual void set_x(T value) noexcept {
            x = value;
        }
        virtual void set_y(T value) noexcept {
            y = value;
        }
        virtual T get_x() const noexcept {
            return x;
        }
        virtual T get_y() const noexcept {
            return y;
        }
        virtual std::pair<T, T> to_pair() {
            return std::pair<T, T>{this->x, this->y};
        }
        virtual long double magnitude() const noexcept override {
            return std::sqrt(static_cast<long double>(x)*static_cast<long double>(x) + // best to avoid call to
                             static_cast<long double>(y)*static_cast<long double>(y)); // std::pow() where possible
        }
        virtual vector2D<long double> unit_vector() {
            return *this / this->magnitude();
        }
        virtual vector2D<T> &rotate(long double angle_in_rad = PI) {

            return *this;
        }
        const T &operator[](unsigned char index) const noexcept override {
            if (index > 1) {
                throw std::invalid_argument("Only the indices '0' and '1' are possible.");
            }
            if (index == 0) {
                return x;
            }
            return y;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(vector2D<U> other) {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(vector2D<U> other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(vector2D<U> other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(vector2D<U> other) {
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator%=(vector2D<U> other) {
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (isIntegral<T> && isIntegral<U>)
        vector2D<T> &operator<<=(vector2D<U> other) {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (isIntegral<T> && isIntegral<U>)
        vector2D<T> &operator>>=(vector2D<U> other) {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (isIntegral<T> && isIntegral<U>)
        vector2D<T> &operator|=(vector2D<U> other) {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (isIntegral<T> && isIntegral<U>)
        vector2D<T> &operator&=(vector2D<U> other) {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (isIntegral<T> && isIntegral<U>)
        vector2D<T> &operator^=(vector2D<U> other) {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        vector2D<T> operator~() requires isIntegral<T> {
            return vector2D<T>(~this->x, ~this->y);
        }
        virtual vector2D<T> operator++(int) noexcept {
            vector2D<T> retvec{this->x, this->y};
            ++this->x;
            ++this->y;
            return retvec;
        }
        virtual vector2D<T> &operator++() noexcept override {
            ++this->x;
            ++this->y;
            return *this;
        }
        virtual vector2D<T> operator--(int) noexcept {
            vector2D<T> retvec{this->x, this->y};
            --this->x;
            --this->y;
            return retvec;
        }
        virtual vector2D<T> &operator--() noexcept override {
            --this->x;
            --this->y;
            return *this;
        }
        virtual vector2D<T> &operator=(const vector<T> &other) noexcept override {
            const vector2D<T> oth = (const vector2D<T> &) (other);
            this->x = oth.x;
            this->y = oth.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator=(const vector2D<U> &other) {
            this->x = other.x;
            this->y = other.y;
        }
        template <isFund U>
        friend std::ostream &operator<<(std::ostream &out, vector2D<U> vec);
        template <isFund U, isFund V>
        friend auto operator+(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isFund U, isFund V>
        friend auto operator-(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isFund U, isFund V>
        friend auto operator*(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x * vec2.x)>;
        template <isFund U, isFund V>
        friend auto operator/(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x / vec2.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator%(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x % vec2.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator<<(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x << vec2.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator>>(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x >> vec2.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator|(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x | vec2.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator&(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x & vec2.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator^(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x ^ vec2.x)>;
        template <isFund U, isFund V>
        friend auto operator+(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x + value)>;
        template <isFund U, isFund V>
        friend auto operator-(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x - value)>;
        template <isFund U, isFund V>
        friend auto operator*(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x * value)>;
        template <isFund U, isFund V>
        friend auto operator/(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x / value)>;
        template <isIntegral U, isIntegral V>
        friend auto operator%(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x % value)>;
        template <isIntegral U, isIntegral V>
        friend auto operator<<(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x << value)>;
        template <isIntegral U, isIntegral V>
        friend auto operator>>(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x >> value)>;
        template <isIntegral U, isIntegral V>
        friend auto operator|(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x | value)>;
        template <isIntegral U, isIntegral V>
        friend auto operator&(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x & value)>;
        template <isFund U, isFund V>
        friend auto operator+(V value, vector2D<U> vec1) -> vector2D<decltype(value + vec1.x)>;
        template <isFund U, isFund V>
        friend auto operator-(V value, vector2D<U> vec1) -> vector2D<decltype(value - vec1.x)>;
        template <isFund U, isFund V>
        friend auto operator*(V value, vector2D<U> vec1) -> vector2D<decltype(value * vec1.x)>;
        template <isFund U, isFund V>
        friend auto operator/(V value, vector2D<U> vec1) -> vector2D<decltype(value / vec1.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator%(V value, vector2D<U> vec1) -> vector2D<decltype(value % vec1.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator<<(V value, vector2D<U> vec1) -> vector2D<decltype(value << vec1.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator>>(V value, vector2D<U> vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator|(V value, vector2D<U> vec1) -> vector2D<decltype(value | vec1.x)>;
        template <isIntegral U, isIntegral V>
        friend auto operator&(V value, vector2D<U> vec1) -> vector2D<decltype(value & vec1.x)>;
        template <isFund U>
        friend class vector2D;
    };
    // template <isFund T, isFund U> // requires (std::is_base_of<vector<typename vec1::vector_type>, vec1>::value && std::is_base_of<vector<typename vec2::vector_type>, vec2>::value)
    // class vector_op_wrapper {
    //     vector2D<T> &v1;
    //     vector2D<U> &v2;
    //     typedef decltype(std::declval<T>() + std::declval<U>()) (*vec_op_func)(T, U);
    // public:
    //     vector_op_wrapper() = delete;
    //     vector_op_wrapper(const vector2D<T> &vec1, const vector2D<U> &vec2) : v1{vec1}, v2{vec2} {}
    //     template <typename funcType>
    //     auto operator()(funcType func) {
    //         if constexpr (sizeof(T) == sizeof(U)) {
    //             if constexpr (std::is_floating_point<T>::value)
    //                 return vector2D<T>(func(v1.get_x(), v2.get_x()), func(v1.get_y(), v2.get_y()));
    //             return vector2D<U>(func(v1.get_x(), v2.get_x()), func(v1.get_y(), v2.get_y()));
    //         }
    //         if constexpr (sizeof(U) > sizeof(T)) {
    //             if constexpr (std::is_integral<U>::value && std::is_floating_point<T>::value)
    //                 return vector2D<T>(func(v1.get_x(), v2.get_x()), func(v1.get_y(), v2.get_y()));
    //             return vector2D<U>(func(v1.get_x(), v2.get_x()), func(v1.get_y(), v2.get_y()));
    //         }
    //         if constexpr (std::is_integral<T>::value && std::is_floating_point<U>::value)
    //             return vector2D<U>(func(v1.get_x(), v2.get_x()), func(v1.get_y(), v2.get_y()));
    //         return vector2D<T>(func(v1.get_x(), v2.get_x()), func(v1.get_y(), v2.get_y()));
    //     }
    // };
    template <isFund U>
    std::ostream &operator<<(std::ostream &out, vector2D<U> vec) {
        return out << +vec.x << "i + " << +vec.y << "j"; // '+' for always printing out numerical value (even for chars)
    }
    template <isFund U, isFund V>
    auto operator+(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x + vec2.x)> {
        return vector2D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y);
    }
    template <isFund U, isFund V>
    auto operator-(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x - vec2.x)> {
        return vector2D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y);
    }
    template <isFund U, isFund V>
    auto operator*(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x * vec2.x)> {
        return vector2D<decltype(vec1.x * vec2.x)>(vec1.x * vec2.x, vec1.y * vec2.y);
    }
    template <isFund U, isFund V>
    auto operator/(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x / vec2.x)> {
        return vector2D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator%(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x % vec2.x)> {
        return vector2D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator<<(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x << vec2.x)> {
        return vector2D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator>>(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x >> vec2.x)> {
        return vector2D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator|(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x | vec2.x)> {
        return vector2D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator&(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x & vec2.x)> {
        return vector2D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator^(vector2D<U> vec1, vector2D<V> vec2) -> vector2D<decltype(vec1.x ^ vec2.x)> {
        return vector2D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y);
    }
    template <isFund U, isFund V>
    auto operator+(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x + value)> {
        return vector2D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value);
    }
    template <isFund U, isFund V>
    auto operator-(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x - value)> {
        return vector2D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value);
    }
    template <isFund U, isFund V>
    auto operator*(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x * value)> {
        return vector2D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value);
    }
    template <isFund U, isFund V>
    auto operator/(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x / value)> {
        return vector2D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value);
    }
    template <isIntegral U, isIntegral V>
    auto operator%(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x % value)> {
        return vector2D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value);
    }
    template <isIntegral U, isIntegral V>
    auto operator<<(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x << value)> {
        return vector2D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value);
    }
    template <isIntegral U, isIntegral V>
    auto operator>>(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x >> value)> {
        return vector2D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value);
    }
    template <isIntegral U, isIntegral V>
    auto operator|(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x | value)> {
        return vector2D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value);
    }
    template <isIntegral U, isIntegral V>
    auto operator&(vector2D<U> vec1, V value) -> vector2D<decltype(vec1.x & value)> {
        return vector2D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value);
    }
    template <isFund U, isFund V>
    auto operator+(V value, vector2D<U> vec1) -> vector2D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y); // but would be an extra func. call
    }
    template <isFund U, isFund V>
    auto operator-(V value, vector2D<U> vec1) -> vector2D<decltype(value - vec1.x)> {
        return vector2D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y);
    }
    template <isFund U, isFund V>
    auto operator*(V value, vector2D<U> vec1) -> vector2D<decltype(value * vec1.x)> {
        return vector2D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y);
    }
    template <isFund U, isFund V>
    auto operator/(V value, vector2D<U> vec1) -> vector2D<decltype(value / vec1.x)> {
        return vector2D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator%(V value, vector2D<U> vec1) -> vector2D<decltype(value % vec1.x)> {
        return vector2D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y);
    }
    template <isIntegral U, isIntegral V> // again, no physical sense, just here for completeness
    auto operator<<(V value, vector2D<U> vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <isIntegral U, isIntegral V> // same as above
    auto operator>>(V value, vector2D<U> vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator|(V value, vector2D<U> vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <isIntegral U, isIntegral V>
    auto operator&(V value, vector2D<U> vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <isFund T> String vector<T>::repr = "0i + 0j + 0k";
}

#endif