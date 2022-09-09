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

template <typename U, typename V>
concept bitwiseOperands = requires (U val1, V val2) {
    {val1 << 1}; // must have all of these operators overloaded
    {val1 >> 1};
    {val1 <<= val2};
    {val1 >>= val2};
    {val1 | val2};
    {val1 & val2};
    {val1 ^ val2};
    {~val1};
    {val1 |= val2};
    {val1 &= val2};
    {val1 ^= val2};
    {val2 <<= val1};
    {val2 >>= val1};
    {val2 | val1};
    {val2 & val1};
    {val2 ^ val1};
    {~val2};
    {val2 |= val1};
    {val2 &= val1};
    {val2 ^= val1};
};

template <typename T>
concept isIntegralNumWrapper = isNumWrapper<T> && requires (T val, T val2, size_t l, long double f) {
    {val%val2} -> isConvertible<T>; // modulo operator
    {val%l} -> isConvertible<T>;
    {l%val} -> isConvertible<T>;
    {val %= val2} -> isConvertible<T>;
    {val %= l} -> isConvertible<T>;
    {val%f} -> isConvertible<T>;
    {f%val} -> isConvertible<T>;
    {val %= f} -> isConvertible<T>;
};

namespace gtd {
    class division_by_zero : public std::invalid_argument {
    public:
        division_by_zero() : std::invalid_argument("Cannot divide by zero.") {}
        division_by_zero(const char *message) : std::invalid_argument(message) {}
    };
    template <isNumWrapper T>
    class vector {
    public:
        vector() = default;
        virtual String str(unsigned char f_p_dec_places = 5) const = 0;
        virtual long double magnitude() const noexcept = 0;
        virtual T &operator[](unsigned char index) = 0;
        virtual const T &operator[](unsigned char index) const = 0;
        virtual vector<T> &operator++() noexcept = 0; //can only declare reference-returning func. for an abstract class
        virtual vector<T> &operator--() noexcept = 0;
        virtual vector<T> &operator=(const vector<T> &other) noexcept = 0; // copies other
        template <isNumWrapper U>
        friend class vector;
    };
    template <isNumWrapper T>
    class vector3D;
    template <isNumWrapper T>
    class vector2D : public vector<T> {
    protected:
        T x{0};
        T y{0};
    public:
        vector2D() = default;
        template <typename U> requires isConvertible<U, T>
        vector2D(const vector2D<U> &other) noexcept : x{other.x}, y{other.y} {}
        template <typename U> requires isConvertible<U, T>
        vector2D(vector2D<U> &&other) noexcept : x{other.x}, y{other.y} {}
        vector2D(T &&x_component, T &&y_component) noexcept : x{x_component}, y{y_component} {}
        vector2D(T &x_component, T &y_component) noexcept : x{x_component}, y{y_component} {}
        String str(unsigned char f_p_dec_places = 5) const override {
            String s;
            s.append_back(x, f_p_dec_places);
            if (y < T{0}) {
                s.append_back("i - ").append_back(-y, f_p_dec_places).append_back("j");
            }
            else {
                s.append_back("i + ").append_back(y, f_p_dec_places).append_back("j");
            }
            return s;
        }
        inline virtual vector2D<T> &set_x(T value) noexcept {
            x = value;
            return *this;
        }
        inline virtual vector2D<T> &set_y(T value) noexcept {
            y = value;
            return *this;
        }
        inline T get_x() const noexcept {
            return x;
        }
        inline T get_y() const noexcept {
            return y;
        }
        inline std::pair<T, T> to_pair() {
            return std::pair<T, T>{this->x, this->y};
        }
        virtual long double magnitude() const noexcept override {
            return std::sqrt(static_cast<long double>(x)*static_cast<long double>(x) + // best to avoid call to
                             static_cast<long double>(y)*static_cast<long double>(y)); // std::pow() where possible
        }
        virtual vector2D<long double> unit_vector() const noexcept {
            if (x == T{0} && y == T{0}) { // no place for a division_by_zero exception if x and y are zero themselves
                return {};
            }
            return *this / this->magnitude();
        }
        virtual vector2D<T> &rotate(long double angle_in_rad = PI) {
            this->apply(matrix<long double>::get_2D_rotate_matrix(PI));
            return *this;
        }
        virtual vector2D<T> &apply(const matrix<T> &transform) {
            if (transform.mat.size() != 2 || transform.mat[0].size() != 2) {
                throw invalid_matrix_format("Only 2x2 matrices can be applied to a vector2D object.");
            }
            x = transform.mat[0][0]*x + transform.mat[0][1]*y;
            y = transform.mat[1][0]*x + transform.mat[1][1]*y;
            return *this;
        }
        inline T &operator[](unsigned char index) noexcept override {
            if (index > 1) {
                throw std::invalid_argument("Only the indices '0' and '1' are possible.");
            }
            if (index == 0) {
                return x;
            }
            return y;
        }
        const T &operator[](unsigned char index) const noexcept override {
            return this->operator[](index);
        }
        virtual vector2D<T> &operator+=(const vector2D<T> &other) { // cannot requires-constrain virtual functions
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        virtual vector2D<T> &operator-=(const vector2D<T> &other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        virtual vector2D<T> &operator*=(const vector2D<T> &other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        virtual vector2D<T> &operator/=(const vector2D<T> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        virtual vector2D<T> &operator+=(const vector2D<T> &&other) {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        virtual vector2D<T> &operator-=(const vector2D<T> &&other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        virtual vector2D<T> &operator*=(const vector2D<T> &&other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        virtual vector2D<T> &operator/=(const vector2D<T> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &other) {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector2D<T> &operator%=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator<<=(const vector2D<U> &other) {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &other) {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &other) {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &other) {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &other) {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &&other) {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &&other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &&other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector2D<T> &operator%=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator<<=(const vector2D<U> &&other) {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &&other) {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &&other) {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &&other) {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &&other) {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
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
            if (&other == this) {
                return *this;
            }
            try {
                const vector2D<T> &oth = dynamic_cast<const vector2D<T>&>(other);
                this->x = oth.x;
                this->y = oth.y;
            } catch (const std::bad_cast &bc) {} // no action taken in case of std::bad_cast - vector obj. is unchanged
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator=(const vector2D<U> &other) {
            this->x = other.x;
            this->y = other.y;
            return *this;
        }
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const vector2D<U> &vec);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x * vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x & value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value % vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator<<(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator>>(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator|(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator&(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector2D<V> &v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x * vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x & value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value % vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator<<(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator>>(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator|(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator&(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector2D<V> &&v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x * vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x & value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value % vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator<<(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator>>(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator|(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator&(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &&m, const vector2D<V> &v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const vector2D<U> &&vec);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x * vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x & value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value % vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator<<(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator>>(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator|(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator&(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &&m, const vector2D<V> &&v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U>
        friend class vector2D;
        template <isNumWrapper U>
        friend class vector3D;
    };
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const vector2D<U> &vec) {
        if (vec.y >= 0) {
            return out << +vec.x << "i + " << +vec.y << "j"; // '+' for always printing out numerical value
        }
        return out << +vec.x << "i - " << +(-vec.y) << "j";
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)> {
        return vector2D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)> {
        return vector2D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x * vec2.x)> {
        return vector2D<decltype(vec1.x * vec2.x)>(vec1.x * vec2.x, vec1.y * vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x << vec2.x)> {
        return vector2D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x >> vec2.x)> {
        return vector2D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x | vec2.x)> {
        return vector2D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x & vec2.x)> {
        return vector2D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x ^ vec2.x)> {
        return vector2D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x + value)> {
        return vector2D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x - value)> {
        return vector2D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x * value)> {
        return vector2D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x << value)> {
        return vector2D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x >> value)> {
        return vector2D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x | value)> {
        return vector2D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x & value)> {
        return vector2D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value - vec1.x)> {
        return vector2D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value * vec1.x)> {
        return vector2D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // again, no physical sense, just here for completeness
    auto operator<<(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // same as above
    auto operator>>(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator|(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator&(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m, const vector2D<V> &v)
    -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
            requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)> {
        return vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>(v).apply(m);
    }
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const vector2D<U> &&vec) {
        if (vec.y >= 0) {
            return out << +vec.x << "i + " << +vec.y << "j"; // '+' for always printing out numerical value
        }
        return out << +vec.x << "i - " << +(-vec.y) << "j";
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x + vec2.x)> {
        return vector2D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x - vec2.x)> {
        return vector2D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x * vec2.x)> {
        return vector2D<decltype(vec1.x * vec2.x)>(vec1.x * vec2.x, vec1.y * vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x << vec2.x)> {
        return vector2D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x >> vec2.x)> {
        return vector2D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x | vec2.x)> {
        return vector2D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x & vec2.x)> {
        return vector2D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x ^ vec2.x)> {
        return vector2D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x + value)> {
        return vector2D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x - value)> {
        return vector2D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x * value)> {
        return vector2D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x << value)> {
        return vector2D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x >> value)> {
        return vector2D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x | value)> {
        return vector2D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x & value)> {
        return vector2D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value - vec1.x)> {
        return vector2D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value * vec1.x)> {
        return vector2D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // again, no physical sense, just here for completeness
    auto operator<<(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // same as above
    auto operator>>(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator|(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator&(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m, const vector2D<V> &&v)
    -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)> {
        return vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>(v).apply(m);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)> {
        return vector2D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)> {
        return vector2D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x * vec2.x)> {
        return vector2D<decltype(vec1.x * vec2.x)>(vec1.x * vec2.x, vec1.y * vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x << vec2.x)> {
        return vector2D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x >> vec2.x)> {
        return vector2D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x | vec2.x)> {
        return vector2D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x & vec2.x)> {
        return vector2D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x ^ vec2.x)> {
        return vector2D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x + value)> {
        return vector2D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x - value)> {
        return vector2D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x * value)> {
        return vector2D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x << value)> {
        return vector2D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x >> value)> {
        return vector2D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x | value)> {
        return vector2D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x & value)> {
        return vector2D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value + vec1.x)> {
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value - vec1.x)> {
        return vector2D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value * vec1.x)> {
        return vector2D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // again, no physical sense, just here for completeness
    auto operator<<(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // same as above
    auto operator>>(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator|(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator&(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &&m, const vector2D<V> &v)
    -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)> {
        return vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>(v).apply(m);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x + vec2.x)> {
        return vector2D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x - vec2.x)> {
        return vector2D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x * vec2.x)> {
        return vector2D<decltype(vec1.x * vec2.x)>(vec1.x * vec2.x, vec1.y * vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x << vec2.x)> {
        return vector2D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x >> vec2.x)> {
        return vector2D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x | vec2.x)> {
        return vector2D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x & vec2.x)> {
        return vector2D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x ^ vec2.x)> {
        return vector2D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x + value)> {
        return vector2D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x - value)> {
        return vector2D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x * value)> {
        return vector2D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x << value)> {
        return vector2D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x >> value)> {
        return vector2D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x | value)> {
        return vector2D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x & value)> {
        return vector2D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value - vec1.x)> {
        return vector2D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value * vec1.x)> {
        return vector2D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0}) {
            throw division_by_zero();
        }
        return vector2D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // again, no physical sense, just here for completeness
    auto operator<<(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V> // same as above
    auto operator>>(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator|(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator&(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &&m, const vector2D<V> &&v)
    -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)> {
        return vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>(v).apply(m);
    }
    template <isNumWrapper T>
    class vector3D : public vector2D<T> {
    private:
        T z{0};
    public:
        vector3D() noexcept = default;
        template <typename U> requires isConvertible<U, T>
        vector3D(const vector2D<U> &other) noexcept : vector2D<T>(other) {}
        template <typename U> requires isConvertible<U, T>
        vector3D(vector2D<U> &&other) noexcept : vector2D<T>(other) {}
        template <typename U> requires isConvertible<U, T>
        vector3D(const vector3D<U> &other) noexcept : vector2D<T>(other), z{other.z} {}
        template <typename U> requires isConvertible<U, T>
        vector3D(vector3D<U> &&other) noexcept : vector2D<T>(other), z{other.z} {}
        vector3D(T x_component, T y_component, T z_component) noexcept : vector2D<T>{x_component, y_component},
        z{z_component} {}
        String str(unsigned char f_p_dec_places = 5) const override {
            String s;
            s.append_back(this->x, f_p_dec_places);
            if (vector2D<T>::y < 0) {
                s.append_back("i - ").append_back(-this->y, f_p_dec_places).append_back("j");
            }
            else {
                s.append_back("i + ").append_back(this->y, f_p_dec_places).append_back("j");
            }
            if (z < 0) {
                s.append_back(" - ").append_back(-this->z, f_p_dec_places).append_back("k");
            }
            else {
                s.append_back(" + ").append_back(this->z, f_p_dec_places).append_back("k");
            }
            return s;
        }
        inline vector3D<T> &set_x(T value) noexcept override {
            this->x = value;
            return *this;
        }
        inline vector3D<T> &set_y(T value) noexcept override {
            this->y = value;
            return *this;
        }
        inline vector3D<T> &set_z(T value) noexcept override {
            this->z = value;
            return *this;
        }
        inline T get_z() {
            return this->z;
        }
        long double magnitude() const noexcept override {
            return std::sqrt(static_cast<long double>(this->x)*static_cast<long double>(this->x) +
                             static_cast<long double>(this->y)*static_cast<long double>(this->y) +
                             static_cast<long double>(this->z)*static_cast<long double>(this->z));
        }
        vector2D<long double> unit_vector() const noexcept {
            return *this / this->magnitude();
        }
        vector2D<T> &rotate(long double angle_in_rad = PI) override {
            this->apply(matrix<long double>::get_2D_rotate_matrix(PI));
            return *this;
        }
        virtual vector3D<T> &apply(const matrix<T> &transform) {
            if (transform.mat.size() != 3 || transform.mat[0].size() != 3) {
                throw invalid_matrix_format("Only 3x3 matrices can be applied to a vector3D object.");
            }
            this->x = transform.mat[0][0]*this->x + transform.mat[0][1]*this->y + transform.mat[0][2]*this->z;
            this->y = transform.mat[1][0]*this->x + transform.mat[1][1]*this->y + transform.mat[1][2]*this->z;
            this->z = transform.mat[2][0]*this->x + transform.mat[2][1]*this->y + transform.mat[2][2]*this->z;
            return *this;
        }
        T &operator[](unsigned char index) noexcept override {
            if (index > 2) {
                throw std::invalid_argument("Only the indices '0', '1' and '2' are possible.");
            }
            if (index == 0) {
                return this->x;
            }
            else if (index == 1) {
                return this->y;
            }
            return this->z;
        }
        const T &operator[](unsigned char index) const noexcept override {
            return this->operator[](index);
        }
        vector3D<T> &operator+=(const vector2D<T> &other) override {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        vector3D<T> &operator-=(const vector2D<T> &other) override {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        vector3D<T> &operator*=(const vector2D<T> &other) override {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0}; // since this is a scalar product (z component in 2D vector is implicitly zero)
            return *this;
        }
        vector3D<T> &operator/=(const vector2D<T> &other) noexcept override {
            this->x /= other.x; // again, vector division makes no physical sense
            this->y /= other.y; // these non-physical methods are only here for programming convenience
            return *this;
        }
        vector3D<T> &operator+=(const vector2D<T> &&other) noexcept override {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        vector3D<T> &operator-=(const vector2D<T> &&other) noexcept override {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        vector3D<T> &operator*=(const vector2D<T> &&other) noexcept override {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0};
            return *this;
        }
        vector3D<T> &operator/=(const vector2D<T> &&other) noexcept override {
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &other) {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(const vector2D<U> &other) {
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector2D<T> &operator%=(const vector2D<U> &other) {
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator<<=(const vector2D<U> &other) {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &other) {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &other) {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &other) {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &other) {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &&other) {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &&other) {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &&other) {
            this->x *= other.x;
            this->y *= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator/=(const vector2D<U> &&other) {
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector2D<T> &operator%=(const vector2D<U> &&other) {
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator<<=(const vector2D<U> &&other) {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &&other) {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &&other) {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &&other) {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &&other) {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        vector3D<T> operator++(int) noexcept {
            vector3D<T> retvec{this->x, this->y, this->z};
            ++this->x;
            ++this->y;
            ++this->z;
            return retvec;
        }
        vector3D<T> &operator++() noexcept override {
            ++this->x;
            ++this->y;
            ++this->z;
            return *this;
        }
        vector3D<T> operator--(int) noexcept {
            vector3D<T> retvec{this->x, this->y, this->z};
            --this->x;
            --this->y;
            --this->z;
            return retvec;
        }
        vector3D<T> &operator--() noexcept override {
            --this->x;
            --this->y;
            --this->z;
            return *this;
        }
        vector3D<T> &operator=(const vector<T> &other) noexcept override {
            if (&other == this) {
                return *this;
            }
            try {
                const vector3D<T> &oth = dynamic_cast<const vector3D<T>&>(other);
                this->x = oth.x;
                this->y = oth.y;
                this->z = oth.z;
            } catch (const std::bad_cast &bc) { // if cast to vector3D fails, try to vector2D
                try {
                    const vector2D<T> &oth_v2 = dynamic_cast<const vector2D<T>&>(other);
                    this->x = oth_v2.x;
                    this->y = oth_v2.y;
                } catch (std::bad_cast &bc2) {}
            }
            return *this;
        }
        vector3D<T> &operator=(const vector3D<T> &other) noexcept {
            if (&other == this) {
                return *this;
            }
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator=(const vector3D<U> &other) {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            return *this;
        }
        template <isNumWrapper U>
        friend class vector2D;
        template <isNumWrapper U>
        friend class vector3D;
    };
}
#endif