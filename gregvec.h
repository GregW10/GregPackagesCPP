//
// Created by mario on 01/09/2022.
//

#ifndef GREGVEC_H
#define GREGVEC_H

#include <iostream>
#include <type_traits>
#include <utility>
#include "gregstring.h"
#include <cmath>
#include <sstream>
#include <any>

template <typename T>
concept isFund = std::is_fundamental<T>::value;

template <typename From, typename To>
concept isConvertible = std::is_convertible<From, To>::value;

template <typename T>
concept isIntegral = std::is_integral<T>::value;

namespace gtd {
    constexpr long double PI = 3.14159265358979323846264338327950288419716939937510582097494459230;
    template <typename... pack> // emulates 'zip()' in Python (I really like what I've done here if I may say so myself)
    std::vector<std::tuple<pack...>> zip(const std::vector<pack>&... vectors) {
        if constexpr (sizeof...(pack) == 0) {
            return std::vector<std::tuple<unsigned char>>();
        }
        std::vector<std::tuple<pack...>> zipped;
        std::vector<size_t> sizes = {(vectors.size())...};
        size_t size = sizes[0];
        if (!std::all_of(sizes.begin(), sizes.end(), [&size](size_t s){return s == size;})) {
            return zipped;
        }
        for (size_t i = 0; i < size; ++i) {;
            zipped.push_back(std::tuple<pack...>(vectors[i]...));
        }
        return zipped;
    }
    template <isFund T, size_t rows1, size_t rows2, size_t columns1, size_t columns2>
    std::vector<std::vector<T>> matmul(T mat1[rows1][columns1], T mat2[rows2][columns2]) {
        std::vector<std::vector<T>> ret_matrix;
        if constexpr (columns1 != rows2 || rows1 == 0 || columns1 == 0 || rows2 == 0 || columns2 == 0) {
            return ret_matrix;
        }
        std::vector<T> sub;
        size_t total;
        for (size_t i = 0; i < rows1; ++i) {
            sub.clear();
            for (size_t j = 0; j < columns2; ++j) {
                total = 0;
                for (size_t k = 0; k < columns1; ++k) {
                    total += mat1[i][k]*mat2[k][j];
                }
                sub.push_back(total); // single element in final matrix
            }
            ret_matrix.push_back(sub); // row in final matrix
        }
        return ret_matrix;
    }
    template <isFund T, isFund U>
    auto matmul(const std::vector<std::vector<T>> &mat1,
                const std::vector<std::vector<U>> &mat2) ->
                std::vector<std::vector<decltype(std::declval<T>()*std::declval<U>())>> {
        std::vector<std::vector<T>> ret_matrix;
        if (mat1.empty() || mat2.empty() || mat1[0].empty() || mat2[0].empty()) {
            return ret_matrix;
        }
        const size_t rows1 = mat1.size();
        const size_t columns1 = mat1[0].size();
        const size_t rows2 = mat2.size();
        const size_t columns2 = mat2[0].size();
        for (const auto &[e1, e2] : zip(mat1, mat2)) { // probing for valid matrix
            if (e1.size() != columns1 || e1.size() != rows2 || e2.size() != columns2) {
                return ret_matrix;
            }
        }
        std::vector<T> sub;
        size_t total;
        for (size_t i = 0; i < rows1; ++i) {
            sub.clear();
            for (size_t j = 0; j < columns2; ++j) {
                total = 0;
                for (size_t k = 0; k < columns1; ++k) {
                    total += mat1[i][k]*mat2[k][j];
                }
                sub.push_back(total);
            }
            ret_matrix.push_back(sub);
        }
        return ret_matrix;
    }
    template <typename T>
    std::ostream &print_matrix_vector(const std::vector<std::vector<T>> &matrix) {
        if (matrix.empty() || matrix[0].empty()) {
            return std::cout;
        }
        std::ostringstream out;
        size_t columns = matrix[0].size();
        size_t count;
        for (const std::vector<T> &col : matrix) {
            out << "|";
            count = 0;
            for (const T &elem : col) {
                out << elem << ", ";
                ++count;
            }
            if (count != columns) {
                return std::cout;
            }
            out.seekp(-2, std::ios_base::end);
            out << "|\n";
        }
        return std::cout << out.str();
    }
    template <isFund T>
    class matrix {
    private:
        std::vector<std::vector<T>> mat;
        void check_validity() {
            size_t column_num;
            if (mat.empty() || (column_num = mat[0].size()) == 0) {
                throw std::invalid_argument("Matrix cannot have zero size.");
            }
            for (const std::vector<T> &row : mat) {
                if (row.size() != column_num) {
                    throw std::invalid_argument("Matrix must have equal number of elements in every row.");
                }
            }
            if (mat.size() == column_num) {
                is_square = true;
            }
        }
        bool is_square = false;
    public:
        matrix() { // constructs a 2x2 identity matrix
            mat.push_back(std::vector<T>{1, 0});
            mat.push_back(std::vector<T>{0, 1});
            is_square = true;
        }
        matrix(size_t rows, size_t columns) {
            if (rows == columns) {
                is_square = true;
            }
            std::vector<T> row;
            for (size_t i = 0; i < rows; ++i) {
                row.clear();
                for (size_t j = 0; j < columns; ++j) {
                    row.push_back(0);
                }
                mat.push_back(row);
            }
        }
        matrix(std::vector<std::vector<T>> &&matrix_elements) :
        mat(std::forward<std::vector<std::vector<T>>>(matrix_elements)) {
            check_validity();
        }
        matrix(const std::vector<std::vector<T>> &matrix_elements) : mat(matrix_elements) {
            check_validity();
        }
        // template <typename U, size_t rows, size_t columns> requires isConvertible<U, T>
        // matrix(const U (&array)[rows][columns]) {
        //     if (rows == columns) {
        //         is_square = true;
        //     }
        //     std::vector<T> row;
        //     for (const U &row : array) {
        //         row.clear();
        //         for (const U &elem : row) {
        //             row.push_back(elem);
        //         }
        //         mat.push_back(row);
        //     }
        // }
        matrix<T> &set_dimensions(size_t num_rows, size_t num_cols);
        size_t num_rows() {
            return mat.size();
        }
        size_t num_cols() {
            return mat[0].size();
        }
        T &at(size_t row_num, size_t col_num) {
            if (row_num > mat.size() || col_num > mat[0].size()) {
                throw std::invalid_argument("Index out of range.");
            }
            return mat[row_num][col_num];
        }
        const T &at(size_t row_num, size_t col_num) const {
            if (row_num > mat.size() || col_num > mat[0].size()) {
                throw std::invalid_argument("Index out of range.");
            }
            return mat[row_num][col_num];
        }
        matrix<T> &make_zero() noexcept {
            for (std::vector<T> &row : mat) {
                for (T &elem : row) {
                    elem = 0;
                }
            }
            return *this;
        }
        matrix<T> &make_identity() noexcept {
            if (!is_square) {
                return *this;
            }
            size_t outer_count = 0;
            size_t inner_count;
            for (std::vector<T> &row : mat) {
                inner_count = 0;
                for (T &elem : row) {
                    elem = inner_count++ == outer_count ? 1 : 0;
                }
                ++outer_count;
            }
            return *this;
        }
        auto begin() const noexcept {
            return mat.cbegin();
        }
        auto end() const noexcept {
            return mat.cend();
        }
        auto begin() noexcept {
            return mat.begin();
        }
        auto end() noexcept {
            return mat.end();
        }
        auto cbegin() const noexcept {
            return mat.cbegin();
        }
        auto cend() const noexcept {
            return mat.cend();
        }
        std::vector<T> operator[](size_t index) const {
            if (index >= mat.size()) {
                throw std::invalid_argument("Index out of range.");
            }
            return {mat[index]};
        }
        template <isFund U>
        friend std::ostream &operator<<(std::ostream &out, const matrix<U> &mat);
        template <isFund U, isFund V>
        friend auto operator*(const matrix<U> &m1, const matrix<V> &m2) ->
        matrix<decltype(std::declval<U>()*std::declval<V>() + std::declval<U>()*std::declval<V>())>;
        template <isFund U>
        friend class matrix;
    };
    template <isFund U>
    std::ostream &operator<<(std::ostream &out, const matrix<U> &mat) {
        for (const std::vector<U> &row : mat) {
            out << '|';
            for (const U &elem : row) {
                out << elem << ' ';
            }
            out << "\b|\n";
        }
        return out;
    }
    template <isFund U, isFund V>
    auto operator*(const matrix<U> &m1, const matrix<V> &m2) ->
    matrix<decltype(std::declval<U>()*std::declval<V>() + std::declval<U>()*std::declval<V>())> {
        using T = decltype(std::declval<U>()*std::declval<V>() + std::declval<U>()*std::declval<V>());
        std::vector<std::vector<T>> r;
        const size_t rows1 = m1.mat.size();
        const size_t columns1 = m1.mat[0].size();
        const size_t rows2 = m2.mat.size();
        const size_t columns2 = m2.mat[0].size();
        for (const auto &[e1, e2] : zip(m1.mat, m2.mat)) { // probing for valid matrix
            std::cout << e1.size() << " " << columns1 << " " << e2.size() << " " << rows2 << " " << columns2 << std::endl;
            if (e1.size() != columns1 || e1.size() != rows2 || e2.size() != columns2) {
                throw std::invalid_argument("For matrix multiplication to be possible, the number of columns in the\n"
                                            "first matrix must equal the number of rows in the second.");
            }
        }
        std::vector<T> sub;
        size_t total;
        for (size_t i = 0; i < rows1; ++i) {
            sub.clear();
            for (size_t j = 0; j < columns2; ++j) {
                total = 0;
                for (size_t k = 0; k < columns1; ++k) {
                    total += m1.at(i, k)*m2.at(k, j); // faster than calling operator[]
                }
                sub.push_back(total);
            }
            r.push_back(sub);
        }
        return matrix<T>(std::move(r));
    }
    template <isFund U, isFund V>
    auto operator+(const matrix<U> &m1, const matrix<V> &m2) ->
    matrix<decltype(std::declval<U>() + std::declval<V>())> {
        using T = decltype(std::declval<U>() + std::declval<V>());
        if (m1.mat.size() != m2.mat.size() || m1.mat[0].size() != m2.mat[0].size()) {
            throw std::invalid_argument("Matrix dimensions must be equal for matrix addition to be possible.");
        }
    }
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