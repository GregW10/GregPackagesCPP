//
//
//

#ifndef GREGDATE_H
#define GREGDATE_H

#include "gregstring.h"
#include <regex>

namespace gtd {
    class Time {

    };
    class Date {
    private:
        String date;
        inline static void check_date(const char *dat) {
            std::regex rgx(R"(^\d\d/\d\d/\d{4}$)");
            if (!std::regex_match(dat, rgx)) {
                throw std::invalid_argument("The date provided is not in the dd/mm/yyyy or mm/dd/yyyy format.");
            }
        }
    public:
        Date() : date("01/01/2000") {}
        Date(const char *dat, bool mmddyyyy = false) {
            check_date(dat);
            date = dat;
        }
        friend std::ostream &operator<<(std::ostream &os, const Date &dat);
        friend Date &operator-(const Date &dat1, const Date &dat2); // yields Time object
    };
    std::ostream &operator<<(std::ostream &os, const Date &dat) {
        os << dat.date;
        return os;
    }
}
#endif