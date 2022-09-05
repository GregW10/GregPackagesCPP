//
// Created by mario on 01/09/2022.
//

#include <iostream>
#include <type_traits>
#include <utility>
#include "gregstring.h"

#ifndef GREGVEC_H
#define GREGVEC_H

template <typename T>
concept isFund = std::is_fundamental<T>::value;

template <typename From, typename To>
concept isConvertible = std::is_convertible<From, To>::value;

namespace gtd {
    template <isFund T>
    class vector {
    public:
        vector(vector<T> &other) noexcept = 0; // copy ctor - no move ctor since class only contains primitive types
        virtual String to_str() const noexcept = 0;
        virtual const char *c_str() const noexcept = 0;
        virtual T &operator[](unsigned char index) noexcept = 0;
        virtual vector &operator++() noexcept = 0; // can only declare reference-returning func. for an abstract class
        virtual vector &operator--() noexcept = 0;
        virtual vector &operator=(vector &&other) = 0; // moves other into current instance
        virtual vector &operator=(vector &other) = 0; // copies other
        virtual T magnitude() = 0;
        // friend std::ostream &operator<<(std::ostream &out, vector<T> &&vec);
        // friend vector operator+(vector &&vec1, vector &&vec2);
        // friend vector operator-(vector &&vec1, vector &&vec2);
        // friend vector operator*(vector &&vec1, vector &&vec2);
        // friend vector operator/(vector &&vec1, vector &&vec2);
    };
    template <isFund T>
    class vector2D : public vector<T> {
    protected:
        T x = 0;
        T y = 0;
        char c_str_buff[70]{0}; // maximum num. chars required to represent a vector3D string
    public:
        unsigned char floating_point_decimal_places = 15; // used for determining number of dec. places in str repr.
        vector2D() = default;
        template <isConvertible<T> U>
        explicit vector2D(vector<U> &other) noexcept : x{other.x}, y{other.y} {}
        vector2D(T &&x_component, T &&y_component) noexcept : x{x_component}, y{y_component} {}
        String to_str() {
            String result;
        }
        T &operator[](unsigned char index) override {
            if (index > 1) {
                throw std::invalid_argument("Only the indices '0' and '1' are possible.");
            }
            if (index == 0) {
                return x;
            }
            return y;
        }
    };
}

#endif