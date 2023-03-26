#ifdef _WIN32
#error "Can't compile on Windows."
#endif

#include <iostream>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <set>
#include <cstdint>
#include <cinttypes>

bool is_numeric_c(char c) {
    return c >= 48 && c <= 57;
}
uint64_t to_ull(const std::string& str) {
    size_t index = 0;
    for (const char& c : str) {
        if (is_numeric_c(c)) {
            break;
        }
        ++index;
    }
    const char *data = str.data() + index;
    uint64_t val = 0;
    while (*data && is_numeric_c(*data)) {
        val *= 10;
        val += *data++ - 48;
    }
    return val;
}
uint64_t to_ull(const char *str) {
    if (!str || !*str) {
        fprintf(stderr, "Non-string!\n");
        abort();
    }
    uint64_t val = 0;
    size_t count = 0;
    while (*str) {
        if (!is_numeric_c(*str)) {
            fprintf(stderr, "Found non-numeric character.\n");
            abort();
        }
        val *= 10;
        val += *str++ - 48;
        ++count;
    }
    if (!count) {
        fprintf(stderr, "Error: could not convert any part of the string to number.\n");
        abort();
    }
    return val;
}
int main(int argc, char **argv) {
    if (argc > 4 || argc < 2) {
        fprintf(stderr, "Error: invalid number of arguments. Usage: cpbmp <freq> [opt: prefix] "
                        "[opt: wait_time]\n");
        return 1;
    }
    uint64_t freq = to_ull(*(argv + 1));
    uint64_t wait_time = argc == 4 ? to_ull(*(argv + 3)) : 600;
    if (!freq) {
        fprintf(stderr, "Error: cannot have zero frequency.\n");
        abort();
    }
    char dir_path[4096];
    if (getcwd(dir_path, 4096) == nullptr) {
        fprintf(stderr, "Error: could not obtain current working directory.\n");
        return 1;
    }
    DIR *dir;
    if ((dir = opendir(dir_path)) == nullptr) {
        fprintf(stderr, "Error: could not open current working directory stream.\n");
        return 1;
    }
    struct dirent* entry;
    std::string path;
    std::set<uint64_t> nums;
    uint64_t num;
    struct stat buff{};
    while (true) {
        while ((entry = readdir(dir)) != nullptr) {
            path = entry->d_name;
            if (path.ends_with(".bmp")) {
                if (((num = to_ull(path)) - 1) % freq || nums.contains(num)) {
                    continue;
                }
                if (stat(path.c_str(), &buff) == -1 || !S_ISREG(buff.st_mode)) {
                    fprintf(stderr, "Error: .bmp file is not valid.\n");
                    abort();
                }
// #ifdef __linux__
  //              while (buff.st_atimespec.)
//#else
                if (time(nullptr) - buff.st_atime < 60) {
                    sleep(300); // to ensure .bmp has finished being written
                }
                std::cout << "Copying \"" << path << "\" with gsutil." << std::endl;
//#endif
                nums.insert(num);
                if (argc == 2) {
                    std::system(("gsutil cp " + path + " gs://nbod_bucket/").c_str());
                } else {
                    std::string command = "gsutil cp ";
                    command.append(path);
                    command += " gs://nbod_bucket/";
                    command.append(*(argv + 2));
                    command.append(path);
                    std::system(command.c_str());
                }
            }
        }
        sleep((unsigned int) wait_time);
        dir = opendir(dir_path);
    }
    return 0;
}
