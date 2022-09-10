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
        vector() noexcept = default;
        vector(const vector<T> &other) noexcept {}
        vector(const vector<T> &&other) noexcept {}
        template <typename U> requires isNumWrapper<U>
        vector(const vector<U> &other) noexcept {}
        template <typename U> requires isNumWrapper<U>
        vector(const vector<U> &&other) noexcept {}
        virtual String str(unsigned char f_p_dec_places = 5) const = 0;
        virtual long double magnitude() const noexcept = 0;
        virtual T &operator[](unsigned char index) = 0;
        virtual const T &operator[](unsigned char index) const = 0;
        virtual vector<T> &operator++() noexcept = 0; //can only declare reference-returning func. for an abstract class
        virtual vector<T> &operator--() noexcept = 0;
        virtual vector<T> &operator=(const vector<T> &other) noexcept = 0; // copies other
        virtual vector<T> &operator=(const vector<T> &&other) noexcept = 0; // copies other
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
        vector2D() noexcept = default;
        vector2D(const vector2D<T> &other) : x{x}, y{y} {}
        vector2D(const vector2D<T> &&other) : x{x}, y{y} {}
        template <typename U> requires isConvertible<U, T>
        vector2D(const vector2D<U> &other) noexcept : x{other.x}, y{other.y} {}
        template <typename U> requires isConvertible<U, T>
        vector2D(const vector2D<U> &&other) noexcept : x{other.x}, y{other.y} {}
        vector2D(const T &x_component, const T &y_component) noexcept : x{x_component}, y{y_component} {}
        vector2D(const T &&x_component, const T &&y_component) noexcept : x{x_component}, y{y_component} {}
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
        virtual vector2D<T> &operator+=(const vector2D<T> &other)noexcept {//cannot requires-constrain virtual functions
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        virtual vector2D<T> &operator-=(const vector2D<T> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        virtual vector2D<T> &operator*=(const vector2D<T> &other) noexcept {
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
        vector2D<T> &operator%=(const vector2D<T> &other) requires isIntegralNumWrapper<T> { // cannot be virtual
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        virtual vector2D<T> &operator+=(const vector2D<T> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        virtual vector2D<T> &operator-=(const vector2D<T> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        virtual vector2D<T> &operator*=(const vector2D<T> &&other) noexcept {
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
        vector2D<T> &operator%=(const vector2D<T> &&other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &other) noexcept {
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
        vector2D<T> &operator<<=(const vector2D<U> &other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator+=(const vector2D<U> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator-=(const vector2D<U> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector2D<T> &operator*=(const vector2D<U> &&other) noexcept {
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
        vector2D<T> &operator<<=(const vector2D<U> &&other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator>>=(const vector2D<U> &&other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator|=(const vector2D<U> &&other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator&=(const vector2D<U> &&other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector2D<T> &operator^=(const vector2D<U> &&other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }
        vector2D<T> operator++(int) noexcept { // cannot make virtual due to differing return types
            vector2D<T> retvec{this->x, this->y};
            ++this->x;
            ++this->y;
            return retvec;
        }
        vector2D<T> &operator++() noexcept override {
            ++this->x;
            ++this->y;
            return *this;
        }
        vector2D<T> operator--(int) noexcept { // cannot make virtual due to differing return types
            vector2D<T> retvec{this->x, this->y};
            --this->x;
            --this->y;
            return retvec;
        }
        vector2D<T> &operator--() noexcept override {
            --this->x;
            --this->y;
            return *this;
        }
        vector2D<T> operator~() requires requires (T val) {~val;}{
            return vector2D<T>(~this->x, ~this->y);
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
        vector2D<T> &operator=(const vector2D<U> &other) noexcept {
            this->x = other.x;
            this->y = other.y;
            return *this;
        }
        vector2D<T> &operator=(const vector<T> &&other) noexcept override {
            return *this = other;
        }
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const vector2D<U> &vec);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x ^ value)>;
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector2D<V> &v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x);
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x ^ value)>;
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector2D<V> &&v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector2D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x ^ value)>;
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value ^ vec1.x)>;
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
        friend auto operator*(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x);
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x ^ value)>;
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
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &&m, const vector2D<V> &&v)
        -> vector2D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector3D<V> &&vec2);
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
    auto operator*(const vector2D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const V &value) -> vector2D<decltype(vec1.x ^ value)> {
        return vector2D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value);
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &value, const vector2D<U> &vec1) -> vector2D<decltype(value ^ vec1.x)> {
        return vector2D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y);
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
    auto operator*(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const V &&value) -> vector2D<decltype(vec1.x ^ value)> {
        return vector2D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value);
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &value, const vector2D<U> &&vec1) -> vector2D<decltype(value ^ vec1.x)> {
        return vector2D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y);
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
    auto operator*(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &&vec1, const V &value) -> vector2D<decltype(vec1.x ^ value)> {
        return vector2D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value);
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &&value, const vector2D<U> &vec1) -> vector2D<decltype(value ^ vec1.x)> {
        return vector2D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y);
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
    auto operator*(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &&vec1, const V &&value) -> vector2D<decltype(vec1.x ^ value)> {
        return vector2D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value + vec1.x)> {
        return vector2D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y);
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
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value << vec1.x)> {
        return vector2D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value >> vec1.x)> {
        return vector2D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value | vec1.x)> {
        return vector2D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value & vec1.x)> {
        return vector2D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &&value, const vector2D<U> &&vec1) -> vector2D<decltype(value ^ vec1.x)> {
        return vector2D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y);
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
        vector3D(const vector2D<T> &other) noexcept : vector2D<T>(other) {}
        vector3D(const vector2D<T> &&other) noexcept : vector2D<T>(other) {}
        vector3D(const vector3D<T> &other) noexcept : vector2D<T>(other), z{other.z} {}
        vector3D(const vector3D<T> &&other) noexcept : vector2D<T>(other), z{other.z} {}
        template <typename U> requires isConvertible<U, T>
        vector3D(const vector2D<U> &other) noexcept : vector2D<T>(other) {}
        template <typename U> requires isConvertible<U, T>
        vector3D(const vector2D<U> &&other) noexcept : vector2D<T>(other) {}
        template <typename U> requires isConvertible<U, T>
        vector3D(const vector3D<U> &other) noexcept : vector2D<T>(other), z{other.z} {}
        template <typename U> requires isConvertible<U, T>
        vector3D(const vector3D<U> &&other) noexcept : vector2D<T>(other), z{other.z} {}
        vector3D(const T &x_component, const T &y_component, const T &z_component) noexcept :
        vector2D<T>{x_component, y_component}, z{z_component} {}
        vector3D(const T &&x_component, const T &&y_component, const T &&z_component) noexcept :
        vector2D<T>{x_component, y_component}, z{z_component} {}
        vector3D(const T &x_component, const T &y_component) noexcept : vector2D<T>{x_component, y_component} {}
        vector3D(const T &&x_component, const T &&y_component) noexcept : vector2D<T>{x_component, y_component} {}
        String str(unsigned char f_p_dec_places = 5) const override {
            String s;
            s.append_back(this->x, f_p_dec_places);
            if (this->y < 0) {
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
        inline vector3D<T> &set_z(T value) noexcept {
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
        vector3D<T> &operator+=(const vector2D<T> &other) noexcept override {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        vector3D<T> &operator-=(const vector2D<T> &other) noexcept override {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        vector3D<T> &operator*=(const vector2D<T> &other) noexcept override {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0}; // given this is the scalar product, mul. of a 3D vec. by a 2D vec. yields a 2D vec (or a 3D
            return *this; // vec. with a zero z-component)
        }
        vector3D<T> &operator/=(const vector2D<T> &other) override {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x; // division between vectors makes no physical sense anyway, so there is no reason to
            this->y /= other.y; // 'divide' by the zero z-component of the 2D vector
            return *this; // these 'non-physical' methods are only here for programming convenience... it should be
        }                 // remembered that they do not represent valid mathematical operations
        vector3D<T> &operator%=(const vector2D<T> &other) requires isIntegralNumWrapper<T> { // hides the parent method,
            if (this->x == T{0} || this->y == T{0}) { // but nothing can be done, since virtual functions cannot have
                throw division_by_zero();             // requires clauses
            }
            this->x %= other.x;
            this->y %= other.y;
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
        vector3D<T> &operator/=(const vector2D<T> &&other) override {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        vector3D<T> &operator%=(const vector2D<T> &&other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector2D<U> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector2D<U> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector2D<U> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0};
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector2D<U> &other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector2D<U> &other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y; // z << 0 would just equal z, so not included
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector2D<U> &other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y; // z >> 0 would just equal z, so not included
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector2D<U> &other) noexcept {
            this->x |= other.x;
            this->y |= other.y; // z | 0 would equal z, so not included
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector2D<U> &other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z = T{0}; // bitwise AND with zero will always equal zero
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector2D<U> &other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y; // z ^ 0 would equal zero, so left unchanged
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector2D<U> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector2D<U> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector2D<U> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z = T{0};
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector2D<U> &&other) {
            if (this->x == T{0} || this->y == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector2D<U> &&other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector2D<U> &&other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector2D<U> &&other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector2D<U> &&other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z = T{0};
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector2D<U> &&other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            return *this;
        }





        vector3D<T> &operator+=(const vector3D<T> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        vector3D<T> &operator-=(const vector3D<T> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        vector3D<T> &operator*=(const vector3D<T> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        vector3D<T> &operator/=(const vector3D<T> &other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        vector3D<T> &operator%=(const vector3D<T> &other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        vector3D<T> &operator+=(const vector3D<T> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        vector3D<T> &operator-=(const vector3D<T> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        vector3D<T> &operator*=(const vector3D<T> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        vector3D<T> &operator/=(const vector3D<T> &&other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        vector3D<T> &operator%=(const vector3D<T> &&other) requires isIntegralNumWrapper<T> {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector3D<U> &other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector3D<U> &other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector3D<U> &other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector3D<U> &other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector3D<U> &other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector3D<U> &other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            this->z <<= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector3D<U> &other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            this->z >>= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector3D<U> &other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            this->z |= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector3D<U> &other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z &= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector3D<U> &other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            this->z ^= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator+=(const vector3D<U> &&other) noexcept {
            this->x += other.x;
            this->y += other.y;
            this->z += other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator-=(const vector3D<U> &&other) noexcept {
            this->x -= other.x;
            this->y -= other.y;
            this->z -= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator*=(const vector3D<U> &&other) noexcept {
            this->x *= other.x;
            this->y *= other.y;
            this->z *= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator/=(const vector3D<U> &&other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x /= other.x;
            this->y /= other.y;
            this->z /= other.z;
            return *this;
        }
        template <typename U> requires (isConvertible<U, T> && isIntegralNumWrapper<T>)
        vector3D<T> &operator%=(const vector3D<U> &&other) {
            if (this->x == T{0} || this->y == T{0} || this->z == T{0}) {
                throw division_by_zero();
            }
            this->x %= other.x;
            this->y %= other.y;
            this->z %= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator<<=(const vector3D<U> &&other) noexcept {
            this->x <<= other.x;
            this->y <<= other.y;
            this->z <<= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator>>=(const vector3D<U> &&other) noexcept {
            this->x >>= other.x;
            this->y >>= other.y;
            this->z >>= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator|=(const vector3D<U> &&other) noexcept {
            this->x |= other.x;
            this->y |= other.y;
            this->z |= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator&=(const vector3D<U> &&other) noexcept {
            this->x &= other.x;
            this->y &= other.y;
            this->z &= other.z;
            return *this;
        }
        template <typename U> requires (bitwiseOperands<T, U>)
        vector3D<T> &operator^=(const vector3D<U> &&other) noexcept {
            this->x ^= other.x;
            this->y ^= other.y;
            this->z ^= other.z;
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
                } catch (std::bad_cast &bc2) {} // this would only occur if a user has created a subclass of vector<T>
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
        vector3D<T> &operator=(const vector3D<U> &other) noexcept {
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
            return *this;
        }
        vector3D<T> &operator=(const vector<T> &&other) noexcept override {
            return *this = other;
        }
        vector3D<T> &operator=(const vector3D<T> &&other) noexcept {
            return *this = other;
        }
        template <typename U> requires (isConvertible<U, T>)
        vector3D<T> &operator=(const vector3D<U> &&other) noexcept {
            return *this = other;
        }
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const vector3D<U> &vec);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x & value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x ^ value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value % vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector3D<V> &v)
        -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x & value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x ^ value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value % vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &m, const vector3D<V> &&v)
        -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x & value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x ^ value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value % vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &&m, const vector3D<V> &v)
        -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>;
        template <isNumWrapper U>
        friend std::ostream &operator<<(std::ostream &out, const vector3D<U> &&vec);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x + value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x - value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x * value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x / value)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x % value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x << value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x >> value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x | value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x & value)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x ^ value)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value + vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value - vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value * vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value / vec1.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value % vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value << vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value >> vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value | vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value & vec1.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value ^ vec1.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const matrix<U> &&m, const vector3D<V> &&v)
        -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
        requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator+(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator-(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator*(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x);
        template <isNumWrapper U, isNumWrapper V>
        friend auto operator/(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)>;
        template <isIntegralNumWrapper U, isIntegralNumWrapper V>
        friend auto operator%(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator<<(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator>>(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator|(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator&(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)>;
        template <typename U, typename V> requires bitwiseOperands<U, V>
        friend auto operator^(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V>
        friend auto cross(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)>;
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector2D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector2D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector2D<U> &&vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector3D<V> &vec2);
        template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
        friend long double angle_between(const vector3D<U> &&vec1, const vector3D<V> &&vec2);
        template <isNumWrapper U>
        friend class vector2D;
        template <isNumWrapper U>
        friend class vector3D;
    };
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const vector3D<U> &vec) {
        out << vec.x << "i ";
        if (vec.y < 0) {
            out << "- " << +(-vec.y) << "j ";
        }
        else {
            out << "+ " << +vec.y << "j ";
        }
        if (vec.z < 0) {
            out << "- " << +(-vec.z) << "k ";
        }
        else {
            out << "+ " << +vec.z << "k ";
        }
        return out;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z / vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z % vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z << vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z >> vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z | vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, vec1.z & vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z ^ vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x + value)> {
        return vector3D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value, vec1.z + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x - value)> {
        return vector3D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value, vec1.z - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x * value)> {
        return vector3D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value, vec1.z * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value, vec1.z / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value, vec1.z % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x << value)> {
        return vector3D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value, vec1.z << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x >> value)> {
        return vector3D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value, vec1.z >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x | value)> {
        return vector3D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value, vec1.z | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x & value)> {
        return vector3D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value, vec1.z & value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const V &value) -> vector3D<decltype(vec1.x ^ value)> {
        return vector3D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value, vec1.z ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector3D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y, value + vec1.z); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value - vec1.x)> {
        return vector3D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y, value - vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value * vec1.x)> {
        return vector3D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y, value * vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y, value / vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y, value % vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value << vec1.x)> {
        return vector3D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y, value << vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value >> vec1.x)> {
        return vector3D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y, value >> vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value | vec1.x)> {
        return vector3D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y, value | vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value & vec1.x)> {
        return vector3D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y, value & vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &value, const vector3D<U> &vec1) -> vector3D<decltype(value ^ vec1.x)> {
        return vector3D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y, value ^ vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m, const vector3D<V> &v)
    -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)> {
        return vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>(v).apply(m);
    }
    template <isNumWrapper U>
    std::ostream &operator<<(std::ostream &out, const vector3D<U> &&vec) {
        return out << vec;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z / vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z % vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z << vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z >> vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z | vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, vec1.z & vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z ^ vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x + value)> {
        return vector3D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value, vec1.z + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x - value)> {
        return vector3D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value, vec1.z - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x * value)> {
        return vector3D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value, vec1.z * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value, vec1.z / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value, vec1.z % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x << value)> {
        return vector3D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value, vec1.z << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x >> value)> {
        return vector3D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value, vec1.z >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x | value)> {
        return vector3D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value, vec1.z | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x & value)> {
        return vector3D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value, vec1.z & value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const V &&value) -> vector3D<decltype(vec1.x ^ value)> {
        return vector3D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value, vec1.z ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value + vec1.x)> { // could just return vec op value,
        return vector3D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y, value + vec1.z); // but would be an extra func. call
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value - vec1.x)> {
        return vector3D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y, value - vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value * vec1.x)> {
        return vector3D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y, value * vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y, value / vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y, value % vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value << vec1.x)> {
        return vector3D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y, value << vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value >> vec1.x)> {
        return vector3D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y, value >> vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value | vec1.x)> {
        return vector3D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y, value | vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value & vec1.x)> {
        return vector3D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y, value & vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &value, const vector3D<U> &&vec1) -> vector3D<decltype(value ^ vec1.x)> {
        return vector3D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y, value ^ vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &m, const vector3D<V> &&v)
    -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)> {
        return vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>(v).apply(m);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z / vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z % vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z << vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z >> vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z | vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, vec1.z & vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z ^ vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x + value)> {
        return vector3D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value, vec1.z + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x - value)> {
        return vector3D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value, vec1.z - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x * value)> {
        return vector3D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value, vec1.z * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value, vec1.z / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value, vec1.z % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x << value)> {
        return vector3D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value, vec1.z << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x >> value)> {
        return vector3D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value, vec1.z >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x | value)> {
        return vector3D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value, vec1.z | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x & value)> {
        return vector3D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value, vec1.z & value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &&vec1, const V &value) -> vector3D<decltype(vec1.x ^ value)> {
        return vector3D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value, vec1.z ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value + vec1.x)> {
        return vector3D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y, value + vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value - vec1.x)> {
        return vector3D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y, value - vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value * vec1.x)> {
        return vector3D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y, value * vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y, value / vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y, value % vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value << vec1.x)> {
        return vector3D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y, value << vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value >> vec1.x)> {
        return vector3D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y, value >> vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value | vec1.x)> {
        return vector3D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y, value | vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value & vec1.x)> {
        return vector3D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y, value & vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &&value, const vector3D<U> &vec1) -> vector3D<decltype(value ^ vec1.x)> {
        return vector3D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y, value ^ vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &&m, const vector3D<V> &v)
    -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)> {
        return vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>(v).apply(m);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z / vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0} || vec2.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z % vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z << vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z >> vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z | vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, vec1.z & vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z ^ vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x + value)> {
        return vector3D<decltype(vec1.x + value)>(vec1.x + value, vec1.y + value, vec1.z + value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x - value)> {
        return vector3D<decltype(vec1.x - value)>(vec1.x - value, vec1.y - value, vec1.z - value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x * value)> {
        return vector3D<decltype(vec1.x * value)>(vec1.x * value, vec1.y * value, vec1.z * value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x / value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / value)>(vec1.x / value, vec1.y / value, vec1.z / value);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x % value)> {
        if (value == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % value)>(vec1.x % value, vec1.y % value, vec1.z % value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x << value)> {
        return vector3D<decltype(vec1.x << value)>(vec1.x << value, vec1.y << value, vec1.z << value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x >> value)> {
        return vector3D<decltype(vec1.x >> value)>(vec1.x >> value, vec1.y >> value, vec1.z >> value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x | value)> {
        return vector3D<decltype(vec1.x | value)>(vec1.x | value, vec1.y | value, vec1.z | value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x & value)> {
        return vector3D<decltype(vec1.x & value)>(vec1.x & value, vec1.y & value, vec1.z & value);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &&vec1, const V &&value) -> vector3D<decltype(vec1.x ^ value)> {
        return vector3D<decltype(vec1.x ^ value)>(vec1.x ^ value, vec1.y ^ value, vec1.z ^ value);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value + vec1.x)> {
        return vector3D<decltype(value + vec1.x)>(value + vec1.x, value + vec1.y, value + vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value - vec1.x)> {
        return vector3D<decltype(value - vec1.x)>(value - vec1.x, value - vec1.y, value - vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value * vec1.x)> {
        return vector3D<decltype(value * vec1.x)>(value * vec1.x, value * vec1.y, value * vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value / vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value / vec1.x)>(value / vec1.x, value / vec1.y, value / vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value % vec1.x)> {
        if (vec1.x == V{0} || vec1.y == V{0} || vec1.z == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(value % vec1.x)>(value % vec1.x, value % vec1.y, value % vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value << vec1.x)> {
        return vector3D<decltype(value << vec1.x)>(value << vec1.x, value << vec1.y, value << vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value >> vec1.x)> {
        return vector3D<decltype(value >> vec1.x)>(value >> vec1.x, value >> vec1.y, value >> vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value | vec1.x)> {
        return vector3D<decltype(value | vec1.x)>(value | vec1.x, value | vec1.y, value | vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value & vec1.x)> {
        return vector3D<decltype(value & vec1.x)>(value & vec1.x, value & vec1.y, value & vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const V &&value, const vector3D<U> &&vec1) -> vector3D<decltype(value ^ vec1.x)> {
        return vector3D<decltype(value ^ vec1.x)>(value ^ vec1.x, value ^ vec1.y, value ^ vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const matrix<U> &&m, const vector3D<V> &&v)
    -> vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>
    requires isConvertible<U, decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)> {
        return vector3D<decltype(std::declval<U>()*v.x + std::declval<U>()*v.y + std::declval<U>()*v.z)>(v).apply(m);
    }






    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec1.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec1.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec1.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator+(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x + vec2.x)> {
        return vector3D<decltype(vec1.x + vec2.x)>(vec1.x + vec2.x, vec1.y + vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator-(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x - vec2.x)> {
        return vector3D<decltype(vec1.x - vec2.x)>(vec1.x - vec2.x, vec1.y - vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator*(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> decltype(vec1.x * vec2.x) {
		return vec1.x*vec2.x + vec1.y*vec2.y;
    }
    template <isNumWrapper U, isNumWrapper V>
    auto operator/(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x / vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x / vec2.x)>(vec1.x / vec2.x, vec1.y / vec2.y, vec2.z);
    }
    template <isIntegralNumWrapper U, isIntegralNumWrapper V>
    auto operator%(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x % vec2.x)> {
        if (vec2.x == V{0} || vec2.y == V{0}) {
            throw division_by_zero();
        }
        return vector3D<decltype(vec1.x % vec2.x)>(vec1.x % vec2.x, vec1.y % vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator<<(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x << vec2.x)> {
        return vector3D<decltype(vec1.x << vec2.x)>(vec1.x << vec2.x, vec1.y << vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator>>(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x >> vec2.x)> {
        return vector3D<decltype(vec1.x >> vec2.x)>(vec1.x >> vec2.x, vec1.y >> vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator|(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x | vec2.x)> {
        return vector3D<decltype(vec1.x | vec2.x)>(vec1.x | vec2.x, vec1.y | vec2.y, vec2.z);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator&(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x & vec2.x)> {
        return vector3D<decltype(vec1.x & vec2.x)>(vec1.x & vec2.x, vec1.y & vec2.y, 0);
    }
    template <typename U, typename V> requires bitwiseOperands<U, V>
    auto operator^(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x ^ vec2.x)> {
        return vector3D<decltype(vec1.x ^ vec2.x)>(vec1.x ^ vec2.x, vec1.y ^ vec2.y, vec2.z);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(0, 0, vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(0, 0, vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(0, 0, vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(0, 0, vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(-vec1.z*vec2.y, // vec2.z is zero
                                                   vec1.z*vec2.x,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(-vec1.z*vec2.y,
                                                   vec1.z*vec2.x,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &&vec1, const vector2D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(-vec1.z*vec2.y,
                                                   vec1.z*vec2.x,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &&vec1, const vector2D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(-vec1.z*vec2.y,
                                                   vec1.z*vec2.x,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z, // vec1.z is zero
                                                   -vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z,
                                                   -vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z,
                                                   -vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector2D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z,
                                                   -vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z - vec1.z*vec2.y,
                                                   vec1.z*vec2.x - vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z - vec1.z*vec2.y,
                                                   vec1.z*vec2.x - vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &&vec1, const vector3D<V> &vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z - vec1.z*vec2.y,
                                                   vec1.z*vec2.x - vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V>
    auto cross(const vector3D<U> &&vec1, const vector3D<V> &&vec2) -> vector3D<decltype(vec1.x * vec2.y)> {
        return vector3D<decltype(vec1.x * vec2.y)>(vec1.y*vec2.z - vec1.z*vec2.y,
                                                   vec1.z*vec2.x - vec1.x*vec2.z,
                                                   vec1.x*vec2.y - vec1.y*vec2.x);
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &vec1, const vector2D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    } // have to check for less than -1 or greater than 1 because of f. p. rounding errors
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &vec1, const vector2D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &&vec1, const vector2D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &&vec1, const vector2D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &vec1, const vector2D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &vec1, const vector2D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &&vec1, const vector2D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &&vec1, const vector2D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &vec1, const vector3D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &vec1, const vector3D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &&vec1, const vector3D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector2D<U> &&vec1, const vector3D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &vec1, const vector3D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &vec1, const vector3D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &&vec1, const vector3D<V> &vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
    template <isNumWrapper U, isNumWrapper V> requires isConvertible<U, long double> && isConvertible<V, long double>
    long double angle_between(const vector3D<U> &&vec1, const vector3D<V> &&vec2) {
        long double cosine_of_angle = (vec1*vec2)/(vec1.magnitude()*vec2.magnitude());
        return cosine_of_angle > 1 ? 0 : (cosine_of_angle < -1 ? PI : std::acos(cosine_of_angle));
    }
}
#endif
