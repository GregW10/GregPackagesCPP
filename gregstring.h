#ifndef GREGSTRING_H
#define GREGSTRING_H

//
// Created by Gregor Hartl Watters on 04/04/2022.
//

#include <cstdlib>
#include <string>
#include <iostream>
#include <iterator>
#include <cstddef>
#include <limits>
#include <vector>
#include <sys/stat.h>
#include <functional>

#ifndef _WIN32
#include <unistd.h>
#include <pwd.h>
#else
#include <shlobj.h>
#endif

static int MIN_SIZE = 128;

namespace gtd {
    unsigned char to_upper(unsigned char ch);
    unsigned char to_lower(unsigned char ch);
    size_t strlen_c(const char *str);
    char *memset_c(char *str, char ch, size_t n_chars);
    char *strcpy_c(char *dest, const char *source);
    inline char *setchr_c(char *str, char ch, size_t pos);
    int strcmp_c(const char *str1, const char *str2);
    int strncmp_c(const char *str1, const char *str2, size_t n);
    bool isdigit_c(const char &ch);

    class NullPointerError : public std::exception {
    private:
        char *message;
    public:
        NullPointerError() {
            const char def_msg[] = "A null pointer cannot be passed as a string.";
            size_t length = strlen_c(def_msg);
            message = new char[length + 1];
            memset_c(message, '\0', length + 1);
            strcpy_c(message, def_msg);
        }

        explicit NullPointerError(const char *msg) {
            if (msg == nullptr) {
                throw std::invalid_argument("nullptr passed as message.");
            }
            size_t length = strlen_c(msg);
            message = new char[length + 1];
            memset_c(message, '\0', length + 1);
            strcpy_c(message, msg);
        }

        [[nodiscard]] const char *what() const noexcept override {
            return message;
        }

        ~NullPointerError() override {
            delete[] message;
        }
    };

    class EmptyStringError : public std::exception {
    private:
        static constexpr char default_msg[] = "This operation cannot be performed on en emtpy string.";
        char *msg = nullptr;
    public:
        EmptyStringError() = default;
        explicit EmptyStringError(const char *message) {
            if (message == nullptr) {
                throw NullPointerError();
            }
            size_t length = strlen_c(message);
            msg = new char[length + 1];
            strcpy_c(msg, message);
            setchr_c(msg, 0, length);
        }
        [[nodiscard]] const char *what() const noexcept override {
            if (msg == nullptr) {
                return default_msg;
            }
            return msg;
        }
        ~EmptyStringError() override {
            delete[] msg;
        }
    };

    class String {
    private:
        char *data = nullptr;
        char *string = nullptr;
        size_t length_w_null = 0;
        size_t size = 0;
        size_t space_front = 0;
        size_t space_back = 0;
        bool is_empty = true;
        bool shrunk = false; // only true after calling shrink_to_fit() and before memory is resized again
        bool start_left = false; // only true if push_left is set to true when a String object is constructed

        static inline bool pow_of_2(size_t num) {
            return (num != 0) && ((num & (num - 1)) == 0); // any power of two will be a 1 followed by zeros in binary,
        } // so, e.g., 10000 & 01111 = 00000 (which is zero) - bitwise AND of n and n-1 for n power of 2 is always zero

        void set_size(bool def_double = true) {
            if (shrunk && !pow_of_2(size)) {
                if (size < MIN_SIZE) {
                    size = MIN_SIZE;
                }
                else {
                    unsigned char n = 0;
                    while ((size >> (++n)) != 0); // find position from RHS of the first non-zero bit (n)
                    size = (~size >> (sizeof(size_t)*8 - n)) + 1; // bitwise negation makes all leading zero bits 1,
                } // then this is bit-shifted by n, which yields 2^n - 1, so +1 after
            }
            if (def_double) {
                size *= 2;
            }
            while (size < length_w_null) {
                size *= 2;
            }
        }

        void constructor(const char *str, bool after_empty = false) {
            size = MIN_SIZE;
            length_w_null = strlen_c(str) + 1;
            set_size(false);
            if (after_empty) { // this is to avoid re-assigning the memory if not necessary
                if (size != MIN_SIZE) {
                    delete[] data;
                    data = new char[size];
                }
            }
            else {
                data = new char[size];
            }
            memset_c(data, '\0', size);
            string = data + get_first_pos();
            // setchr_c(string, '\0', length_w_null - 1);
            strcpy_c(string, str);
            space_front = get_first_pos();
            space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front + 1);
            is_empty = false;
        }

        void empty_constructor() {
            size = MIN_SIZE;
            data = new char[size];
            memset_c(data, '\0', size);
            string = nullptr;
            length_w_null = 0;
            if (start_left) {
                space_front = 0;
                space_back = size - length_w_null;
            }
            else {
                space_front = size / 2;
                space_back = size / 2;
            }
            is_empty = true;
        }

        template <typename It>
        void it_constructor(It beg, It ending) {
            length_w_null = ending - beg + 1;
            size = MIN_SIZE;
            set_size(false);
            data = new char[size];
            memset_c(data, '\0', size);
            string = data + get_first_pos();
            char *ptr = string;
            for (; beg != ending; ++beg) {
                *ptr = *beg;
                ++ptr;
            }
            setchr_c(string, '\0', length_w_null - 1);
            space_front = get_first_pos();
            space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front + 1);
            is_empty = false;
        }

        inline static bool is_even(size_t num) {
            return num % 2 == 0;
        }

        inline static bool is_punc(const char &ch) { // reference because it is a private func. only used with l-values
            return (ch > 0 && ch < 48) || (ch > 57 && ch < 65) || (ch > 90 && ch < 97) || ch > 122;
        }

        [[nodiscard]] inline size_t get_first_pos() const noexcept {
            return start_left ? 0 : (is_even(length_w_null) ? size/2 - length_w_null/2 : size/2 - (length_w_null+1)/2);
        }

    public:
        static const size_t nopos = -1;
        class RevIterator;
        class Iterator {
        public:
            using iterator_category = std::bidirectional_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = char;
            using pointer = char *;
            using reference = char &;
            Iterator() : ptr(nullptr) {}
            Iterator(char *pointer) : ptr(pointer) {}
            reference operator*() const {
                return *ptr;
            }
            pointer operator->() const {
                return ptr;
            }
            virtual Iterator &operator++() {
                ++ptr;
                return *this;
            }
            Iterator operator++(int) {
                Iterator tmp = *this;
                ++ptr;
                return tmp;
            }
            virtual Iterator &operator--() {
                --ptr;
                return *this;
            }
            Iterator operator--(int) {
                Iterator tmp = *this;
                --ptr;
                return tmp;
            }
            friend bool operator==(const Iterator &A, const Iterator &B);
            friend bool operator!=(const Iterator &A, const Iterator &B);
            friend bool operator<(const Iterator &A, const Iterator &B);
            friend bool operator>(const Iterator &A, const Iterator &B);
            friend bool operator<=(const Iterator &A, const Iterator &B);
            friend bool operator>=(const Iterator &A, const Iterator &B);
            Iterator &operator=(const Iterator &other) {
                if (*this == other) {
                    return *this;
                }
                this->ptr = other.ptr;
                return *this;
            }
            Iterator &operator=(char *pointer) {
                if (this->ptr == pointer) {
                    return *this;
                }
                this->ptr = pointer;
                return *this;
            }
            Iterator operator+(std::ptrdiff_t offset) {
                return {ptr + offset};
            }

            Iterator operator-(std::ptrdiff_t offset) {
                return {ptr - offset};
            }
            friend std::ostream &operator<<(std::ostream &out, const Iterator &A);
            friend difference_type operator-(const Iterator &A, const Iterator &B);
            friend std::ptrdiff_t operator-(const String::RevIterator &A, const String::RevIterator &B);
            ~Iterator() = default;
        protected:
            char *ptr;
        };
        class RevIterator : public Iterator {
        public:
            RevIterator() = default;
            RevIterator(char *pointer) : Iterator(pointer) {}
            RevIterator &operator++() override {
                --ptr;
                return *this;
            }
            RevIterator operator++(int) {
                RevIterator tmp = *this;
                --ptr;
                return tmp;
            }
            RevIterator &operator--() override {
                ++ptr;
                return *this;
            }
            RevIterator operator--(int) {
                RevIterator tmp = *this;
                ++ptr;
                return tmp;
            }
            friend bool operator<(const RevIterator &A, const RevIterator &B);
            friend bool operator>(const RevIterator &A, const RevIterator &B);
            friend bool operator<=(const RevIterator &A, const RevIterator &B);
            friend bool operator>=(const RevIterator &A, const RevIterator &B);
            RevIterator operator+(std::ptrdiff_t offset) {
                return {ptr - offset};
            }

            RevIterator operator-(std::ptrdiff_t offset) {
                return {ptr + offset};
            }
        };
        class ConstIterator : public Iterator {
        private:
            const char *cptr;
        public:
            ConstIterator() = default;
            ConstIterator(char *pointer) : Iterator(pointer), cptr(pointer){}
            const char *operator->() {
                cptr = ptr;
                return cptr;
            }
            const char &operator*() {
                cptr = ptr;
                return *cptr;
            }
        };
        class ConstRevIterator : public RevIterator {
        private:
            const char *cptr;
        public:
            ConstRevIterator() : cptr(nullptr) {}
            ConstRevIterator(char *pointer) : RevIterator(pointer), cptr(pointer){}
            const char *operator->() {
                cptr = ptr;
                return cptr;
            }
            const char &operator*() {
                cptr = ptr;
                return *cptr;
            }
        };

        String(bool push_left = false) {
            start_left = push_left;
            empty_constructor();
        }

        String(char ch, bool push_left = false) {
            start_left = push_left;
            if (ch == '\0') {
                empty_constructor();
            }
            else {
                const char str[2]{ch, '\0'};
                constructor(str);
            }
        }

        String(const char *str, bool push_left = false) {
            start_left = push_left;
            if (str == nullptr) {
                throw NullPointerError();
            }
            if (strlen_c(str) == 0) {
                empty_constructor();
            } else {
                constructor(str);
            }
        }

        String(const std::string &str, bool push_left = false) {
            start_left = push_left;
            if (str.empty()) {
                empty_constructor();
            } else {
                constructor(str.c_str());
            }
        }

        String(const String &str, bool push_left = false) {
            start_left = push_left;
            if (str.empty()) {
                empty_constructor();
            } else {
                constructor(str.c_str());
            }
        }

        template <typename ForwardIterator>
        String(ForwardIterator beg, ForwardIterator ending, bool push_left = false) {
            start_left = push_left;
            if (beg >= ending) {
                empty_constructor();
                return;
            }
            it_constructor(beg, ending);
        }

        void append_back(const char *str) {
            if (str == nullptr) {
                throw NullPointerError();
            }
            if (is_empty) {
                constructor(str, true);
                return;
            }
            size_t l = strlen_c(str);
            size_t old_l = length_w_null - 1;
            length_w_null += l;
            if (l > space_back) {
                if (shrunk) {
                    set_size(false);
                }
                else {
                    set_size();
                }
                char *old_data = data;
                char *old_str = string;
                data = new char[size];
                memset_c(data, '\0', size);
                string = data + get_first_pos();
                strcpy_c(string, old_str);
                delete[] old_data;
                space_front = get_first_pos();
                space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front+1);
            }
            strcpy_c(string + old_l, str);
            space_back -= l;
        }

        void append_back(const std::string &str) {
            if (str.empty()) {
                return;
            }
            append_back(str.c_str());
        }

        void append_back(const String &str) {
            if (str.is_empty) {
                return;
            }
            append_back(str.c_str());
        }

        void append_front(const char *str) {
            if (str == nullptr) {
                throw NullPointerError();
            }
            if (is_empty) {
                constructor(str, true);
                return;
            }
            size_t l = strlen_c(str);
            length_w_null += l;
            char gone = *string;
            if (l > space_front) {
                if (shrunk) {
                    set_size(false);
                }
                else {
                    set_size();
                }
                char *old_data = data;
                char *old_str = string;
                data = new char[size];
                memset_c(data, '\0', size);
                string = data + get_first_pos();
                strcpy_c(string + l, old_str);
                delete[] old_data;
                strcpy_c(string, str);
                space_front = get_first_pos();
                space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front : space_front+1);
            } else {
                strcpy_c(string - l, str);
                string -= l;
                space_front -= l;
            }
            setchr_c(string, gone, l);
        }

        void print_all_chars() {
            const char *ptr = data;
            const char *end = data + size;
            size_t count = 0;
            for (; ptr < end; ++ptr) {
                if (*ptr == 0) {
                    printf("0");
                    ++count;
                    continue;
                }
                else {
                    printf("%c", *ptr);
                    ++count;
                }
            }
            printf("\nCount: %zu\n", count);
        }

        void append_front(const std::string &str) {
            if (str.empty()) {
                return;
            }
            append_front(str.c_str());
        }

        void append_front(const String &str) {
            if (str.is_empty) {
                return;
            }
            append_front(str.c_str());
        }

        void push_back(char ch) {
            if (ch == 0) {
                return;
            }
            const char str[2] = {ch, '\0'};
            append_back(str);
        }

        void push_front(char ch) {
            if (ch == 0) {
                return;
            }
            const char str[2] = {ch, '\0'};
            append_front(str);
        }

        bool pop_back() noexcept {
            if (is_empty) {
                return false;
            }
            if (length_w_null == 2 || length_w_null == 1) {
                this->clear();
                return true;
            }
            *(string + length_w_null - 2) = '\0';
            --length_w_null;
            ++space_back;
            return true;
        }

        bool pop_front() noexcept {
            if (is_empty) {
                return false;
            }
            if (length_w_null == 2 || length_w_null == 1) {
                this->clear();
                return true;
            }
            *string++ = '\0';
            --length_w_null;
            ++space_front;
            return true;
        }

        void reverse() noexcept {
            if (is_empty || length_w_null == 1 || length_w_null == 2) {
                return;
            }
            char *start = string;
            char *final = start + length_w_null - 2;
            char temp;
            do {
                temp = *start;
                *start = *final;
                *final = temp;
            } while (--final > ++start);
        }

        void fill(char ch) noexcept {
            if (is_empty || ch == '\0') {
                return;
            }
            char *ptr = string + length_w_null - 2;
            while (ptr >= string) {
                *ptr-- = ch;
            }
        }

        void fill(size_t starting_index, size_t end_index, char ch) noexcept { // including end_index
            if (starting_index >= end_index || starting_index < 0 || end_index < 0 || ch == '\0' || is_empty ||
            starting_index >= length_w_null - 1 || end_index > length_w_null - 1) {
                return;
            }
            char *start = string + starting_index;
            const char *final = end_index > length_w_null - 1 ? string + length_w_null - 1 : string + end_index;
            while (start <= final) {
                *start++ = ch;
            }
        }

        void fill_data(char ch) noexcept {
            if (ch == '\0') {
                return;
            }
            char *start = data;
            const char *final = data + size;
            for (; start != final; ++start) {
                *start = ch;
            }
            *(--start) = '\0';
            string = data;
            length_w_null = size;
            space_front = 0;
            space_back = 0;
            is_empty = false;
        }

        template <typename ForwardIterator>
        void assign(ForwardIterator start, ForwardIterator final) {
            if (start >= final) {
                return;
            }
            delete[] data;
            it_constructor(start, final);
        }

        String substr(size_t starting_index, size_t end_index) const {
            if (starting_index < 0 || end_index < 0 || starting_index >= end_index || starting_index > length_w_null -
                2 || is_empty) {
                return {};
            }
            if (end_index > length_w_null - 1) {
                end_index = length_w_null - 1;
            }
            Iterator start = begin() + (long int) starting_index;
            Iterator final = begin() + (long int) end_index;
            return {start, final};
        }

        String substr(char ch) const noexcept {
            if (ch == '\0' || is_empty) {
                return {};
            }
            char *start = string;
            const char *final = string + length_w_null - 1;
            for (; start != final; ++start) {
                if (*start == ch) {
                    return {start};
                }
            }
            return {};
        }

        String substr(const char *str) const noexcept {
            size_t len = strlen_c(str);
            if (str == nullptr || len == 0 || len > length_w_null - 1 || is_empty) {
                return {};
            }
            if (len == length_w_null - 1) {
                if (strcmp_c(string, str) == 0) {
                    return {string};
                }
                return {};
            }
            char *start_m = string;
            char *start;
            char *final_m = string + len;
            const char *final_fixed = string + length_w_null - 1;
            int count = 0;
            while (final_fixed >= final_m) {
                for (start = start_m; start != final_m; ++start) {
                    if (*start != *(str + count)) {
                        break;
                    }
                    ++count;
                }
                if (count == len) {
                    return {start_m};
                }
                count = 0;
                ++start_m;
                ++final_m;
            }
            return {};
        }

        String rsubstr(char ch) const noexcept {
            if (ch == '\0' || is_empty) {
                return {};
            }
            const char *start = string;
            const char *final = start + length_w_null - 1;
            --final;
            for (; start != final; --final) {
                if (*final == ch) {
                    return {string + (final - start)};
                }
            }
            if (*start == ch) {
                return {string};
            }
            return {};
        }

        String rsubstr(const char *str) const noexcept {
            size_t len = strlen_c(str);
            if (str == nullptr || len == 0 || len > length_w_null - 1 || is_empty) {
                return {};
            }
            if (len == length_w_null - 1) {
                if (strcmp_c(string, str) == 0) {
                    return {string};
                }
                return {};
            }
            const char *start_fixed = string;
            const char *start_m = string + length_w_null - 1 - len;
            const char *start;
            const char *final_m = string + length_w_null - 1;
            size_t count = 0;
            while (start_fixed <= start_m) {
                for (start = start_m; start != final_m; ++start) {
                    if (*start != *(str + count)) {
                        break;
                    }
                    ++count;
                }
                if (count == len) {
                    return {start_m};
                }
                count = 0;
                --start_m;
                --final_m;
            }
            return {};
        }

        String &insert(char ch, size_t pos = nopos) {// moves the chars in whichever direction has the least chars to
            if (is_empty) {                         // move, or in the direction which does not cause a resize in memory
                const char str[2]{ch, 0};
                constructor(str, true);
                return *this;
            }
            if (pos > length_w_null - 2) {
                if (space_back != 0 || (space_front == 0 && space_back == 0)) {
                    push_back(ch);
                    return *this;
                }
                pos = length_w_null - 1;
            }
            if (pos == 0) {
                push_front(ch);
                return *this;
            }
            char *to_place = string + pos;
            if (space_front != 0 || space_back != 0) {
                if (pos >= length_w_null / 2 || (pos < length_w_null / 2 && space_front == 0)) {
                    char *ptr = string + length_w_null - 2;
                    char *end = string + length_w_null - 1;
                    while (ptr >= to_place) {
                        *end-- = *ptr--;
                    }
                    --space_back;
                    *++to_place = ch;
                } else {
                    char *beg = string - 1;
                    char *ptr = string;
                    while (ptr <= to_place) {
                        *beg++ = *ptr++;
                    }
                    --string;
                    --space_front;
                    *--to_place = ch;
                }
                ++length_w_null;
                return *this;
            }

        }

        void to_upper() noexcept {
            if (is_empty) {
                return;
            }
            char *ptr = string + length_w_null - 2;
            while (ptr >= string) {
                if (*ptr >= 97 && *ptr <= 122) {
                    *ptr -= 32;
                }
                --ptr;
            }
        }

        void to_lower() noexcept {
            if (is_empty) {
                return;
            }
            char *ptr = string + length_w_null - 2;
            while (ptr >= string) {
                if (*ptr >= 65 && *ptr <= 90) {
                    *ptr += 32;
                }
                --ptr;
            }
        }

        bool isnumeric() const noexcept {
            if (is_empty) {
                return false;
            }
            char *start = string;
            const char *final = start + length_w_null - 1;
            while (start != final) {
                if (!isdigit_c(*start++)) {
                    return false;
                }
            }
            return true;
        }

        bool contains(char ch) const noexcept {
            if (is_empty) {
                return false;
            }
            char *start = string;
            const char *final = start + length_w_null - 1;
            while (start != final) {
                if (*start++ == ch) {
                    return true;
                }
            }
            return false;
        }

        bool contains(const char *str) const noexcept {
            if (is_empty || str == nullptr || *str == 0 || strlen_c(str) > length_w_null - 1) {
                return false;
            }
            size_t length = strlen_c(str);
            if (length == length_w_null - 1) {
                return *this == str;
            }
            const char *sptr;
            const char *ptr = string;
            const char *end_ptr = string + length_w_null - length;
            const char *vptr;
            size_t i;
            size_t count;
            for (; ptr <= end_ptr; ++ptr) {
                for (count = 0, i = 0, sptr = str, vptr = ptr; i < length; ++i, ++sptr, ++vptr) {
                    if (*sptr != *vptr) {
                        break;
                    }
                    ++count;
                }
                if (count == length) {
                    return true;
                }
            }
            return false;
        }

        size_t strip(char ch) noexcept { // returns the number of characters removed
            if (is_empty || !contains(ch) || ch == '\0') {
                return 0;
            }
            char *start = string;
            size_t temp_len = length_w_null - 1;
            size_t count = 0;
            while (*start != '\0') {
                if (*start == ch) {
                    for (int i = 0; i < temp_len; ++i) {
                        *(start + i) = *(start + i + 1);
                    }
                    --length_w_null;
                    --space_back;
                    start = string;
                    temp_len = length_w_null - 1;
                    ++count;
                    continue;
                }
                ++start;
                --temp_len;
            }
            if (length_w_null < 2) {
                this->clear();
            }
            return count;
        }

        size_t strip(const char *characters) noexcept {
            if (is_empty) {
                return 0;
            }
            size_t len_char = strlen_c(characters);
            size_t not_count = 0;
            for (const char &ch : *this) {
                if (!contains(ch)) {
                    ++not_count;
                }
            }
            if (len_char == not_count) {
                return 0;
            }
            size_t count = 0;
            for (size_t i = 0; i < len_char; ++i) {
                count += strip(*characters++);
            }
            return count;
        }

        size_t remove(char ch) noexcept {
            if (is_empty || !contains(ch) || ch == '\0') {
                return nopos;
            }
            char *ptr = string;
            size_t pos = 0;
            while (true) {
                if (*ptr == ch) {
                    for (; *ptr != '\0'; ++ptr) {
                        *ptr = *(ptr + 1);
                    }
                    --space_back;
                    --length_w_null;
                    if (length_w_null == 1) {
                        this->clear();
                    }
                    return pos;
                }
                ++ptr;
                ++pos;
            }
        }

        size_t remove(const char *str) {
            size_t len = strlen_c(str);
            if (is_empty || str == nullptr || *str == '\0' || length_w_null - 1 < len) {
                return nopos;
            }
            if (length_w_null - 1 == len) {
                if (strcmp_c(string, str) == 0) {
                    this->clear();
                    return 0;
                }
                return nopos;
            }
            char *optr = string;
            char *ptr = string;
            char *sub = new char[len + 1];
            char *optr_sub;
            size_t count;
            for (size_t pos = 0; pos < length_w_null - len; ++pos) {
                count = 0;
                optr_sub = optr;
                while (count < len) {
                    *sub = *optr_sub;
                    ++count; ++sub; ++optr_sub;
                }
                *sub = 0;
                sub -= len;
                if (strcmp_c(str, sub) == 0) {
                    for (size_t i = 0; i < len; ++i) {
                        ptr = optr;
                        for (; *ptr; ++ptr) {
                            *ptr = *(ptr + 1);
                        }
                        --space_back;
                        --length_w_null;
                    }
                    if (length_w_null == 1) {
                        this->clear();
                    }
                    return pos;
                }
                ++optr;
                ++ptr;
            }
            return nopos;
        }

        size_t r_remove(char ch) noexcept {
            if (is_empty || !contains(ch) || ch == '\0') {
                return nopos;
            }
            char *ptr = string + length_w_null - 2;
            size_t pos = length_w_null - 2;
            while (true) {
                if (*ptr == ch) {
                    for (; *ptr != '\0'; ++ptr) {
                        *ptr = *(ptr + 1);
                    }
                    --space_back;
                    --length_w_null;
                    if (length_w_null == 1) {
                        this->clear();
                    }
                    return pos;
                }
                --ptr;
                --pos;
            }
        }

        std::vector<String> split(char delim = 32) const {
            if (is_empty) {
                return {};
            }
            if (!contains(delim)) {
                return {string};
            }
            std::vector<String> retvec;
            char *ptr = string;
            gtd::String part(true);
            while (*ptr) {
                while (*ptr != delim && *ptr) {
                    part.push_back(*ptr);
                    ++ptr;
                }
                if (!part.empty()) {
                    retvec.push_back(part);
                }
                part.clear();
                ++ptr;
            }
            return retvec;
        }

        bool startswith(const char *beg) const noexcept {
            if (!contains(beg)) {
                return false;
            }
            return strncmp_c(string, beg, strlen_c(beg)) == 0;
        }

        bool endswith(const char *end) const noexcept {
            if (!contains(end)) {
                return false;
            }
            const char *str = string + length_w_null - strlen_c(end) - 1;
            return strcmp_c(str, end) == 0;
        }

        int adopt_text(const char *path) noexcept {
            if (path == nullptr || *path == 0) {
                return -1;
            }
            struct stat buffer{};
            FILE *fp;
            if (stat(path, &buffer) == -1 || S_ISDIR(buffer.st_mode) || buffer.st_size == 0 ||
                (fp = fopen(path, "r")) == nullptr || fgetc(fp) == EOF) {
                return -1;
            }
            if (!is_empty) {
                this->clear();
            }
            fseek(fp, 0, SEEK_SET);
            size_t count = 0;
            while (fgetc(fp) != EOF) { // instead of using struct stat, because st_size would include EOF & others
                ++count;
            }
            fseek(fp, 0, SEEK_SET);
            length_w_null = count + 1;
            set_size(false);
            data = new char[size];
            string = data + get_first_pos();
            space_front = get_first_pos();
            space_back = start_left ? size - length_w_null : (is_even(length_w_null) ? space_front+1 : space_front+2);
            is_empty = false;
            memset_c(data, 0, size);
            fread(string, sizeof(char), count, fp);
            fclose(fp);
            return 0;
        }

        size_t word_count() {
            if (is_empty) {
                return 0;
            }
            size_t word_count = 0;
            size_t word_len = 0;
            bool has_backspace = false;
            bool has_carriage_return = false;
            const char *ptr = string;
            while (*ptr) {
                if (*ptr >= 33 && *ptr <= 126) {
                    word_len = 1;
                    ++word_count;
                    while (*(++ptr) >= 33 && *ptr <= 126) {
                        ++word_len;
                    }
                    if (*ptr == 0) {
                        return word_count;
                    }
                    if (*ptr == '\b' && (*(ptr + 1) >= 33 && *(ptr + 1) <= 126)) {
                        --word_count;
                    }
                }
                ++ptr;
            }
        }

        size_t shift_center() noexcept {
            if (is_empty || space_back == space_front || space_back == space_front + 1) {
                return 0;
            }
            size_t move_by = space_back > space_front ? (is_even(length_w_null) ? (space_back - space_front) / 2 :
                    (space_back - space_front - 1) / 2) : (is_even(length_w_null) ? (space_front - space_back) / 2 :
                            (space_front - space_back + 1) / 2);
            char *first_ptr = string;
            char *fixed_first = first_ptr;
            char *last_ptr = string + length_w_null - 1;
            char *less_ptr = nullptr;
            char *fixed_last = last_ptr;
            if (space_back > space_front) {
                for (size_t i = 0; i < move_by; ++i) {
                    last_ptr = fixed_last;
                    less_ptr = last_ptr - 1;
                    for (; last_ptr > first_ptr; --last_ptr, --less_ptr) {
                        *last_ptr = *less_ptr;
                        *less_ptr = 0;
                    }
                    ++first_ptr;
                    ++fixed_last;
                }
                string += move_by;
                space_front += move_by;
                space_back -= move_by;
                return move_by;
            }
            for (size_t i = 0; i < move_by; ++i) {
                first_ptr = fixed_first;
                less_ptr = first_ptr - 1;
                for (; first_ptr <= last_ptr; ++first_ptr, ++less_ptr) {
                    *less_ptr = *first_ptr;
                }
                --fixed_first;
                --last_ptr;
            }
            string -= move_by;
            space_front -= move_by;
            space_back += move_by;
            return move_by;
        }

        size_t shift_left() noexcept {
            if (is_empty || space_front == 0) {
                return 0;
            }
            char *less_ptr = string - 1;
            char *fixed_ptr = string;
            char *ptr = string;
            const char *end = string + length_w_null;
            while (fixed_ptr > data) {
                while (ptr < end) {
                    *less_ptr = *ptr;
                    ++ptr; ++less_ptr;
                }
                --end;
                --fixed_ptr;
                ptr = fixed_ptr;
                less_ptr = fixed_ptr - 1;
            }
            size_t retval = space_front;
            space_front = 0;
            space_back = size - length_w_null;
            string = data;
            return retval;
        }

        size_t shift_right() noexcept {
            if (is_empty || space_back == 0) {
                return 0;
            }
            char *fixed_beg = string;
            char *ptr = string + length_w_null - 1;
            char *less_ptr = ptr - 1;
            char *fixed_end = ptr;
            const char *end = data + size - 1;
            while (fixed_end < end) {
                while (less_ptr >= fixed_beg) {
                    *ptr = *less_ptr;
                    --ptr; --less_ptr;
                }
                *fixed_beg = 0;
                ++fixed_beg; ++fixed_end;
                ptr = fixed_end; less_ptr = ptr - 1;
            }
            string += space_back;
            size_t retval = space_back;
            space_back = 0;
            space_front = size - length_w_null;
            return retval;
        }

        size_t transform_chars(std::function<void(char&)> func) { // returns the number of characters altered
            if (is_empty) { //                 ^^^^^^ func not passed as reference in case of r-value lambda expressions
                return 0;
            }
            char c;
            size_t transformed_cnt = 0;
            for (char &ch : *this) {
                c = ch;
                func(ch);
                if (ch == 0) {
                    ch = c;
                    continue;
                }
                if (c != ch) {
                    ++transformed_cnt;
                }
            }
            return transformed_cnt;
        }

        size_t transform_words(std::function<const char *(char*)> func, bool check_punctuation = true) { // returns the number of words altered
            if (is_empty || !contains(32)) {
                return 0;
            }
            char *start_ptr = string;
            char *ptr = string;
            char *start_of_word = nullptr;
            char *word = nullptr;
            const char *word_copy = nullptr;
            const char *new_word = nullptr;
            size_t word_len = 0;
            size_t words_altered = 0;
            bool end_punc = false;
            bool tripunc = false;
            bool altered = false;
            size_t num_puncs = 0;
            while (*ptr) {
                if (*ptr == 32) {
                    ++ptr;
                    continue;
                }
                end_punc = false;
                tripunc = false;
                altered = false;
                num_puncs = 0;
                start_ptr = ptr;
                start_of_word = ptr;
                word_len = 0;
                while (*ptr && *ptr != 32) {
                    if (is_punc(*ptr)) {
                        ++num_puncs;
                    }
                    ++word_len;
                    ++ptr;
                }
                if (check_punctuation) {
                    if (num_puncs > 0) {
                        if (num_puncs == 1) {
                            if (is_punc(*(ptr - 1)) && word_len > 1) {
                                end_punc = true;
                                --ptr;
                                --word_len;
                            }
                            else if (is_punc(*start_ptr) && word_len > 1) {
                                ++start_ptr; ++start_of_word;
                                --word_len;
                            }
                        }
                        else if (num_puncs == 2) {
                            if (is_punc(*(ptr - 1)) && is_punc(*start_ptr) && word_len > 2) {
                                end_punc = true;
                                --word_len; --word_len;
                                --ptr; ++start_ptr; ++start_of_word;
                            }
                        }
                        else if (num_puncs == 3) {
                            if (*(ptr - 1) == '.' && (*(ptr - 2) == '\'' || *(ptr - 2) == '"') &&
                                (*start_ptr == '\'' || *start_ptr == '"')) {
                                tripunc = true;
                                --word_len; --word_len; --word_len;
                                --ptr; --ptr; ++start_ptr; ++start_of_word;
                            }
                        }
                    }
                }
                word = new char[word_len + 1];
                word_copy = word;
                for (size_t i = 0; i < word_len; ++i, ++start_of_word) {
                    *(word + i) = *start_of_word;
                }
                *(word + word_len) = 0;
                new_word = func(word);
                if (strcmp_c(new_word, word) == 0) { // same word, so same length
                    if (new_word != word) {
                        delete[] word;
                        return nopos;
                    }
                    for (size_t i = 0; start_ptr < ptr; ++start_ptr, ++i) {
                        if (!altered && *start_ptr != *(word + i)) {
                            altered = true;
                            ++words_altered;
                        }
                        if (*(word + i) == 0) {
                            *(word + i) = 32;
                        }
                        *start_ptr = *(word + i);
                    }
                    if (end_punc) {
                        ++ptr;
                    } else if (tripunc) {
                        ++ptr;
                        ++ptr;
                    }
                }
                else {

                    size_t new_len = strlen_c(new_word);
                }
                delete[] word;
            }
            return words_altered;
        }

        size_t shrink_to_fit() {
            if (is_empty || space_front == 0 || space_back == 0) {
                return 0;
            }
            shrunk = true;
            const char *old_ptr = data;
            data = new char[length_w_null];
            strcpy_c(data, string);
            delete[] old_ptr;
            string = data;
            size_t retval = size - length_w_null;
            size = length_w_null;
            space_front = space_back = 0;
            return retval;
        }

        void space() {
            printf("Space front: %zu, space back: %zu, string - data: %d\n", space_front, space_back, string - data);
        }

        void clear() noexcept {
            delete[] data;
            empty_constructor();
        }
        // erase_chars() shifts the string in whichever direction saves computation
        String &erase_chars(size_t start_pos = 0, size_t end_pos = nopos) noexcept { // erases including end_pos
            if (is_empty || start_pos >= length_w_null - 1 || start_pos > end_pos) {
                return *this;
            }
            if (end_pos > length_w_null - 2) {
                end_pos = length_w_null - 2;
            }
            char *str = string + start_pos;
            char *end = string + end_pos;
            char *end_end = string + length_w_null - 1;
            size_t diff = end_pos - start_pos + 1;
            if (start_pos == 0) { // case for entire beginning-of-string being erased
                while (string <= end) {
                    *string++ = 0;
                }
                if (end_pos == length_w_null - 2) { // case for entire string being erased
                    is_empty = true;
                    string = nullptr;
                    if (start_left) {
                        space_front = 0;
                        space_back = size;
                    } else {
                        space_front = space_back = size / 2;
                    }
                    length_w_null = 0;
                    return *this;
                }
                space_front += diff;
                length_w_null -= diff;
                return *this;
            }
            if (end_pos == length_w_null - 2) { // case for entire end-of-string being erased
                while (end >= str) {
                    *end-- = 0;
                }
                space_back += diff;
                length_w_null -= diff;
                return *this;
            }
            length_w_null -= diff;
            if (str - string >= end_end - end) { // case for end-of-string being shifted left
                ++end;
                while (end <= end_end) {
                    *str++ = *end++;
                }
                while (str <= end_end) {
                    *str++ = 0;
                }
                space_back += diff;
                return *this;
            } // below: case for end-of-string being shifted right
            --str;
            while (str >= string) {
                *end-- = *str--;
            }
            while (end >= string) {
                *end-- = 0;
            }
            string += diff;
            space_front += diff;
            return *this;
        }

        Iterator begin() const noexcept {
            if (is_empty) { return {}; }
            return {string};
        }

        Iterator end() const noexcept {
            if (is_empty) { return {}; }
            return {string + length_w_null - 1};
        }

        RevIterator rbegin() const noexcept {
            if (is_empty) { return {}; }
            return {string + length_w_null - 2};
        }

        RevIterator rend() const noexcept {
            if (is_empty) { return {}; }
            return {string - 1};
        }

        ConstIterator cbegin() const noexcept {
            if (is_empty) { return {}; }
            return {string};
        }

        ConstIterator cend() const noexcept {
            if (is_empty) { return {}; }
            return {string + length_w_null - 1};
        }

        ConstRevIterator crbegin() const noexcept {
            if (is_empty) { return {}; }
            return {string + length_w_null - 2};
        }

        ConstRevIterator crend() const noexcept {
            if (is_empty) { return {}; }
            return {string - 1};
        }

        char &front() const {
            if (is_empty) {
                throw EmptyStringError("front() cannot be called on an empty string.");
            }
            return *string;
        }

        char &back() const {
            if (is_empty) {
                throw EmptyStringError("back() cannot be called on an empty string.");
            }
            return *(string + length_w_null - 2);
        }

        size_t get_size() const noexcept {
            return size;
        }

        size_t get_length() const noexcept {
            return length_w_null == 0 ? 0 : length_w_null - 1;
        }

        bool empty() const noexcept {
            return is_empty;
        }

        const char *const c_str() const noexcept {
            return string;
        }

        std::string str() const noexcept {
            return {string};
        }

        ~String() {
            delete[] data;
        }

        friend std::ostream &operator>>(const String &str, std::ostream &os);

        char &operator[](size_t index) const {
            if (is_empty) {
                throw std::out_of_range("No characters to be accessed in an empty string.");
            }
            if (index > length_w_null - 2) {
                throw std::out_of_range("You are indexing a character that is out of range.");
            }
            return *(string + index);
        }

        bool operator==(const char *str) const noexcept {
            if (str == nullptr) {
                return false;
            }
            size_t length = strlen_c(str);
            if (is_empty) {
                return length == 0;
            }
            if (length != length_w_null - 1) {
                return false;
            }
            const char *copy_string = string;
            while (*str) {
                if (*str != *copy_string) {
                    return false;
                }
                ++str;
                ++copy_string;
            }
            return true;
        }

        bool operator==(const String &str) const {
            return *this == str.c_str();
        }

        bool operator==(const std::string &str) const {
            return *this == str.c_str();
        }

        String operator+(const char *str) {
            String retstr = *this;
            retstr.append_back(str);
            return retstr;
        }

        String operator+(const String &str) {
            String retstr = *this;
            retstr.append_back(str);
            return retstr;
        }

        String operator+(const std::string &str) {
            String retstr = *this;
            retstr.append_back(str);
            return retstr;
        }

        String &operator+=(const char *str) {
            *this = *this + str;
            return *this;
        }

        String &operator+=(const std::string &str) {
            *this = *this + str;
            return *this;
        }

        String &operator+=(const String &str) {
            *this = *this + str;
            return *this;
        }

        String &operator=(const char *str) {
            if (str == nullptr) {
                throw NullPointerError();
            }
            if (*this == str) {
                return *this;
            }
            if (is_empty) {
                constructor(str, true);
                return *this;
            }
            this->clear();
            if (strlen_c(str) == 0) {
                return *this;
            }
            this->constructor(str);
            return *this;
        }

        String &operator=(const std::string &str) {
            this->clear();
            if (str.empty()) {
                return *this;
            }
            return *this = str.c_str();
        }

        String &operator=(const String &str) {
            if (this == &str) {
                return *this;
            }
            return *this = str.str();
        }

        String &operator++() {
            if (is_empty) {
                return *this;
            }
            for (char &ch : *this) {
                if (ch < 127) {
                    ++ch;
                }
            }
            return *this;
        }

        String operator++(int) {
            if (is_empty) {
                return {};
            }
            for (char &ch : *this) {
                if (ch < 127) {
                    ++ch;
                }
            }
            String new_str = *this;
            return new_str;
        }

        String &operator--() {
            if (is_empty) {
                return *this;
            }
            for (char &ch : *this) {
                if (ch > 1) {
                    --ch;
                }
            }
            return *this;
        }

        String operator--(int) {
            if (is_empty) {
                return {};
            }
            for (char &ch : *this) {
                if (ch > 1) {
                    --ch;
                }
            }
            String new_str = *this;
            return new_str;
        }
        friend std::istream &operator>>(std::istream &is, String &str);
        friend std::istream &getline(std::istream &is, String &str, char delim);
    };

    std::istream &operator>>(std::istream &is, String &str) {
        if (!is.good()) {
            return is;
        }
        if (is.peek() == EOF) {
            is.setstate(std::ios::failbit);
            return is;
        }
        bool org_s_left = str.start_left;
        if (!org_s_left) {
            str.start_left = true;
        }
        if (!str.is_empty) {
            str.clear();
        }
        char c;
        while (is.good() && (c = is.get()) != 32 && c != '\n') {
            str.push_back(c);
        }
        is.clear();
        str.start_left = org_s_left;
        return is;
    }

    std::istream &getline(std::istream &is, String &str, char delim = '\n') {
        if (!is.good() || is.peek() == EOF) {
            return is;
        }
        bool org_s_left = str.start_left;
        if (!org_s_left) {
            str.start_left = true;
        }
        if (!str.is_empty) {
            str.clear();
        }
        char ch;
        while (is.good() && (ch = is.get()) != delim) {
            str.push_back(ch);
        }
        str.start_left = org_s_left;
        return is;
    }

    bool operator==(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr == B.ptr;
    }

    bool operator!=(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr != B.ptr;
    }

    bool operator<(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr < B.ptr;
    }

    bool operator>(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr > B.ptr;
    }

    bool operator<=(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr <= B.ptr;
    }

    bool operator>=(const String::Iterator &A, const String::Iterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr >= B.ptr;
    }

    std::ptrdiff_t operator-(const String::Iterator &A, const String::Iterator &B) {
        return A.ptr - B.ptr;
    }

    std::ostream &operator<<(std::ostream &out, const String::Iterator &A) {
        out << static_cast<void *>(A.ptr);
        return out;
    }

    std::ostream &operator>>(const String &str, std::ostream &os) {
        if (str.is_empty) {
            return os;
        }
        os << str.string;
        return os;
    }

    bool operator<(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr > B.ptr;
    }
    bool operator>(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr < B.ptr;
    }
    bool operator<=(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr >= B.ptr;
    }
    bool operator>=(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return false;
        }
        return A.ptr <= B.ptr;
    }

    std::ptrdiff_t operator-(const String::RevIterator &A, const String::RevIterator &B) {
        if (A.ptr == nullptr || B.ptr == nullptr) {
            return 0;
        }
        return B.ptr - A.ptr;
    }

    std::ostream &operator<<(std::ostream &out, const String &str) {
        return str >> out;
    }

    bool operator==(const char *str, const String &string) {
        return string == str;
    }

    bool operator==(const std::string &str, const String &string) {
        return string == str;
    }

    String operator+(const char *str, const String &string) {
        String retstr = string;
        retstr.append_front(str);
        return retstr;
    }

    String operator+(const std::string &str, const String &string) {
        String retstr = string;
        retstr.append_front(str);
        return retstr;
    }

    size_t strlen_c(const char *str) {
        if (str == nullptr) {
            throw NullPointerError();
        }
        size_t length_c = 0;
        while (*str) {
            ++length_c;
            ++str;
        }
        return length_c;
    }

    char *memset_c(char *str, char ch, size_t n_chars) {
        if (str == nullptr) {
            return nullptr;
        }
        for (int i = 0; i < n_chars; i++) {
            *(str + i) = ch;
        }
        return str;
    }

    char *strcpy_c(char *dest, const char *source) {
        if (dest == nullptr || source == nullptr) {
            return nullptr;
        }
        char *ptr = dest;
        while ((*(dest++) = *(source++)));
        return ptr;
    }

    int strcmp_c(const char *str1, const char *str2) {
        if (str1 == nullptr || str2 == nullptr) {
            return -128;
        }
        if (strlen_c(str1) != strlen_c(str2)) {
            return -128;
        }
        while (*str1 && *str2) {
            if (*str1 != *str2) {
                return *str1 - *str2;
            }
            ++str1;
            ++str2;
        }
        return 0;
    }

    int strncmp_c(const char *str1, const char *str2, size_t n) {
        if (str1 == nullptr || str2 == nullptr) {
            return -128;
        }
        size_t count = 0;
        while (*str1 && *str2 && count < n) {
            if (*str1 != *str2) {
                return *str1 - *str2;
            }
            ++str1;
            ++str2;
            ++count;
        }
        return 0;
    }

    inline char *setchr_c(char *str, const char ch, size_t pos) {
        if (str == nullptr) {
            return nullptr;
        }
        *(str + pos) = ch;
        return str;
    }

    inline unsigned char to_upper(unsigned char ch) {
        int diff = 97 - 65;
        if (ch <= 122 && ch >= 97) {
            ch -= diff;
        }
        return ch;
    }

    inline unsigned char to_lower(unsigned char ch) {
        int diff = 97 - 65;
        if (ch <= 90 && ch >= 65) {
            ch += diff;
        }
        return ch;
    }

    std::string to_upper(const char *str) {
        size_t length = strlen_c(str);
        std::string ret_string;
        for (int i = 0; i < length; i++) {
            ret_string.push_back((char) to_upper(*(str + i)));
        }
        return ret_string;
    }

    std::string to_upper(const std::string &str) {
        size_t length = str.length();
        std::string ret_string;
        for (const char &ch: str) {
            ret_string.push_back((char) to_upper(ch));
        }
        return ret_string;
    }

    void string_upper(std::string &str) {
        for (char &ch: str) {
            ch = (char) to_upper(ch);
        }
    }

    void string_upper(char *str) {
        if (str == nullptr) {
            return;
        }
        while (*str) {
            *str = (char) to_upper(*str); ++str;
        }
    }

    void string_lower(std::string &str) {
        for (char &ch: str) {
            ch = (char) to_lower(ch);
        }
    }

    void string_lower(char *str) {
        size_t length = strlen_c(str);
        for (int i = 0; i < length; i++) {
            *(str + i) = (char) to_lower(*(str + i));
        }
    }

    inline bool isdigit_c(const char &ch) { // reference, since this function would never be used for an r-value char
        return ch >= 48 && ch <= 57;
    }

    bool is_numeric(const char *str) {
        if (str == nullptr || *str == 0) {
            return false;
        }
        while (*str) {
            if (!isdigit_c(*str)) {
                return false;
            }
            ++str;
        }
        return true;
    }

    bool is_numeric(const std::string &str) {
        if (str.empty()) {
            return false;
        }
        for (const char &ch: str) {
            if (isdigit_c(ch) == 0) {
                return false;
            }
        }
        return true;
    }

    bool contains(const char *str, char ch) {
        if (str == nullptr || *str == 0 || ch == 0) {
            return false;
        }
        while (*str) {
            if (*str == ch) {
                return true;
            }
            ++str;
        }
        return false;
    }

    size_t count_char(const char *str, char ch = 32) {
        if (str == nullptr || *str == 0) {
            return 0;
        }
        size_t num = 0;
        while (*str) {
            if (*str == ch) {
                ++num;
            }
            ++str;
        }
        return num;
    }

    char **strsplit(const char *str, char delim = 32) {
        if (str == nullptr || *str == 0) {
            return nullptr;
        }
        char **retptr;
        if (delim == 0 || !contains(str, delim)) {
            retptr = new char*[2*sizeof(char *)];
            *(retptr + 1) = nullptr;
            *retptr = new char[strlen_c(str) + 1];
            strcpy_c(*retptr, str);
            return retptr;
        }
        size_t str_len = 0;
        size_t str_count = 0;
        if (*str == delim) {
            ++str;
        }
        const char *beg = str;
        size_t num_strings = count_char(str, delim) + 1;
        retptr = new char*[num_strings + 1];
        while (str_count != num_strings) {
            if (*str == delim || *str == 0) {
                *(retptr + str_count) = new char[str_len + 1];
                for (size_t i = 0; i < str_len; ++i, ++beg) {
                    *(*(retptr + str_count) + i) = *beg;
                }
                *(*(retptr + str_count) + str_len) = 0;
                ++beg;
                ++str_count;
                str_len = 0;
                ++str;
                continue;
            }
            ++str;
            ++str_len;
        }
        *(retptr + num_strings) = nullptr;
        return retptr;
    }

    void free_split_str(const char *const *array) {
        if (array == nullptr || *array == nullptr) {
            return;
        }
        const char *const *ptr = array;
        while (*ptr != nullptr) {
            delete[] *ptr;
            ++ptr;
        }
        delete[] *ptr;
        delete[] array;
    }

    void clear_cin() {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    template<typename PATH>
    PATH get_home_path();

    template<> char *get_home_path<char *>() {
#ifdef _WIN32
        char *home_path_c = (char *) malloc(MAX_PATH);
        HRESULT result = SHGetFolderPathA(nullptr, CSIDL_PROFILE, nullptr, SHGFP_TYPE_CURRENT, home_path_c);
        if (result != S_OK) {
            return nullptr;
        }
        home_path_c = (char *) realloc(home_path_c, strlen_c(home_path_c) + 1);
        return home_path_c;
#else
        struct passwd *pwd;
        uid_t uid = getuid();
        pwd = getpwuid(uid);
        char *home_path_c = pwd->pw_dir;
        char *retval;
        retval = (char *) malloc(sizeof(char) * strlen_c(home_path_c) + 1);
        memset_c(retval, '\0', strlen_c(home_path_c) + 1);
        strcpy_c(retval, home_path_c);
        return retval;
#endif
    }

    template<> std::string get_home_path<std::string>() {
#ifdef _WIN32
        char home_path_c[MAX_PATH];
        HRESULT result = SHGetFolderPathA(nullptr, CSIDL_PROFILE, nullptr, SHGFP_TYPE_CURRENT, home_path_c);
        if (result != S_OK) {
            return {};
        }
#else
        struct passwd *pwd;
        uid_t uid = getuid();
        pwd = getpwuid(uid);
        char *home_path_c = pwd->pw_dir;
#endif
        return {home_path_c};
    }

    template<> String get_home_path<String>() {
#ifdef _WIN32
        char home_path_c[MAX_PATH];
        HRESULT result = SHGetFolderPathA(nullptr, CSIDL_PROFILE, nullptr, SHGFP_TYPE_CURRENT, home_path_c);
        if (result != S_OK) {
            return {};
        }
#else
        struct passwd *pwd;
        uid_t uid = getuid();
        pwd = getpwuid(uid);
        char *home_path_c = pwd->pw_dir;
#endif
        return {home_path_c};
    }

    template<typename whatever>
    void print(whatever w) {
        std::cout << w << std::endl;
    }

    void print(bool w) {
        std::cout << std::boolalpha << w << std::endl;
    }

    template <typename T, typename... types>
    void print(T arg, types... args) {
        print(arg);
        if constexpr (sizeof...(args) > 0) {
            print(args...);
        }
    }

    template<typename T>
    inline void swap(T &A, T &B) {
        T C{A};
        A = B;
        B = C;
    }

    template<class BiDirectionalIterator>
    void reverse(BiDirectionalIterator begin, BiDirectionalIterator end) {
        end--;
        if (begin == end) {
            return;
        }
        while (begin < end) {
            swap(*begin, *end);
            ++begin;
            --end;
        }
    }

    template<class ForwardIterator, typename T>
    void fill(ForwardIterator begin, ForwardIterator end, T value) {
        if (begin == end && begin == end - 1) {
            return;
        }
        while (begin != end) {
            *begin = value;
            ++begin;
        }
    }

    namespace word_transforms {
        typedef const char *(*WT)(char *);
        const char *reverse(char *word) {
            size_t length = strlen_c(word);
            if (word == nullptr || length == 0 || length == 1) {
                return word;
            }
            char *beg = word;
            char *end = word + strlen_c(word) - 1;
            char temp;
            do {
                temp = *beg;
                *beg = *end;
                *end = temp;
            } while (++beg < --end);
            return word;
        }
        const char *cap_first_letter(char *word) {
            if (word == nullptr || *word == 0) {
                return word;
            }
            *word = *word <= 122 && *word >= 97 ? *word - 32 : *word;
            return word;
        }
    }
}
#endif