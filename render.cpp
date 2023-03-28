#include "gregsys.hpp"

gtd::vec3 extract_vec(const std::string& str) {
    long double x, y;
    x = std::stold(str, &index);
    str.erase(0, index);
    y = std::stold(str, &index);
    str.erase(0, index);
    return {x, y, std::stold(str)};
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Invalid number of command-line arguments.\nUsage: ./render <nsys_file> <camera_config>\n";
        return 1;
    }
    std::ifstream in{*(argv + 2)};
    if (!in) {
        std::cerr << "Error opening camera configuration file.\n";
        return 1;
    }
    bool have_cam_info = false;
    std::string str;
    std::getline(in, str);
    if (!std::regex_match(str.c_str(), R"(^\s*(I|i)mage\s*(D|d)imensions\s*=\s*\d{0,10}\s+\d{0,10}\s*$)")) {
        std::cerr << "First line does not match expected format: \"Image dimensions = xxxx xxxx\"";
        return 1;
    }
    gtd::image_dimensions dims{};
    size_t index = str.find('=');
    str.erase(0, index + 1);
    dims.x = std::stoui(str, &index);
    str.erase(0, index);
    dims.y = std::stoui(str);
    do {std::getline(in, str);} while (std::regex_match(str.c_str(), R"(^\s*$)"));
    if (!std::regex_match(str, R"(^\s*(C|c)amera\s*:\s*$)")) {
        std::cerr << "Invalid line: \"" << str << "\" instead of: \"Camera:\"\n";
        return 1;
    }
    std::getline(in, str);
    if (!std::regex_match(str, R"(^\s*(P|p)osition\s*=\s*\d{0,1000}(\.\d{0,1000})?\s+\d{0,1000}(\.\d{0,1000})?\s+\d{0,1000}(\.\d{0,1000})?\s*$)")) {
        std::cerr << "Invalid line: \"" << str << "\" instead of: \"Position = x y z\"\n";
        return 1;
    }
    index = str.find('=');
    str.erase(0, index + 1);
    gtd::vec3 cam_pos = extract_vec(str);
    std::getline(in, str);
    if (!std::regex_match(str, R"(^\s*(D|d)irection\s*=\s*\d{0,1000}(\.\d{0,1000})?\s+\d{0,1000}(\.\d{0,1000})?\s+\d{0,1000}(\.\d{0,1000})?\s*$)")) {
        std::cerr << "Invalid line: \"" << str << "\" instead of: \"Direction = x y z\"\n";
        return 1;
    }
    index = str.find('=');
    str.erase(0, index + 1);
    gtd::vec3 cam_pos = extract_vec(str);
    std::getline(in, str);
    if (!std::regex_match(str, R"(^\s*(R|r)otation\s*(\s*rad\s*)\s*=\s*\d{0,1000}(\.\d{0,1000})?\s*$)")) {
        std::cerr << "Invalid line: \"" << str << "\" instead of: \"Rotation (rad) = x\"\n";
        return 1;
    }
    index = str.find('=');
    str.erase(0, index + 1);
    long double rot = std::stold(str);
    std::getline(in, str);
    if (!std::regex_match(str, R"(^\s*(I|i)mage\s+(D|d)istance\s*=\s*\d{0,1000}(\.\d{0,1000})?\s*$)")) {
        std::cerr << "Invalid line: \"" << str << "\" instead of: \"Rotation (rad) = x\"\n";
        return 1;
    }
    index = str.find('=');
    str.erase(0, index + 1);
    long double im_dist = std::stold(str);
    do {std::getline(in, str);} while (std::regex_match(str.c_str(), R"(^\s*$)"));
    return 0;
}
