#include "simplex/simplex.hpp"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Minimal JSON parser for the specific LP data format:
//   {"name":"...","m":INT,"n":INT,"A":[[...],...],"b":[...],"c":[...],...}
// ---------------------------------------------------------------------------

static std::string read_file(const std::string& path)
{
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open: " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// Find the value of an integer key like "m": 27
static int parse_int_key(const std::string& json, const std::string& key)
{
    // search for "key": digits
    std::string pattern = "\"" + key + "\"";
    std::size_t pos = json.find(pattern);
    if (pos == std::string::npos)
        throw std::runtime_error("Key not found: " + key);
    pos = json.find(':', pos) + 1;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t')) ++pos;
    return std::stoi(json.substr(pos));
}

// Parse a flat JSON number array starting at json[start] which must be '['.
// Returns the numbers and sets end to one past the closing ']'.
static std::vector<double> parse_number_array(const std::string& json,
                                               std::size_t start,
                                               std::size_t& end)
{
    std::vector<double> vals;
    if (json[start] != '[')
        throw std::runtime_error("Expected '[' at position " + std::to_string(start));
    std::size_t pos = start + 1;
    while (pos < json.size()) {
        // skip whitespace / commas
        char c = json[pos];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == ',') {
            ++pos;
            continue;
        }
        if (c == ']') {
            end = pos + 1;
            return vals;
        }
        // parse number
        std::size_t num_end;
        double v = std::stod(json.substr(pos), &num_end);
        vals.push_back(v);
        pos += num_end;
    }
    throw std::runtime_error("Unterminated array");
}

// Parse a 2-D JSON array [[...], [...], ...] starting at json[start].
static std::vector<std::vector<double>> parse_2d_array(const std::string& json,
                                                        std::size_t start,
                                                        std::size_t& end)
{
    std::vector<std::vector<double>> rows;
    if (json[start] != '[')
        throw std::runtime_error("Expected '[' for 2D array");
    std::size_t pos = start + 1;
    while (pos < json.size()) {
        char c = json[pos];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == ',') {
            ++pos;
            continue;
        }
        if (c == ']') {
            end = pos + 1;
            return rows;
        }
        if (c == '[') {
            std::size_t row_end;
            rows.push_back(parse_number_array(json, pos, row_end));
            pos = row_end;
            continue;
        }
        throw std::runtime_error("Unexpected char in 2D array");
    }
    throw std::runtime_error("Unterminated 2D array");
}

// Find the start of the value for a given key (returns index of the value char).
static std::size_t find_value_start(const std::string& json, const std::string& key)
{
    std::string pattern = "\"" + key + "\"";
    std::size_t pos = json.find(pattern);
    if (pos == std::string::npos)
        throw std::runtime_error("Key not found: " + key);
    pos = json.find(':', pos) + 1;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t' ||
                                  json[pos] == '\n' || json[pos] == '\r'))
        ++pos;
    return pos;
}

struct LpData {
    int m, n;
    Eigen::MatrixXd A;
    Eigen::VectorXd b, c;
};

static LpData load_lp(const std::string& path)
{
    std::string json = read_file(path);

    LpData lp;
    lp.m = parse_int_key(json, "m");
    lp.n = parse_int_key(json, "n");

    // Parse A
    std::size_t a_start = find_value_start(json, "A");
    std::size_t a_end;
    auto rows = parse_2d_array(json, a_start, a_end);
    if (static_cast<int>(rows.size()) != lp.m)
        throw std::runtime_error("A row count mismatch");
    lp.A.resize(lp.m, lp.n);
    for (int i = 0; i < lp.m; ++i)
        for (int j = 0; j < lp.n; ++j)
            lp.A(i, j) = rows[i][j];

    // Parse b
    std::size_t b_start = find_value_start(json, "b");
    std::size_t b_end;
    auto bv = parse_number_array(json, b_start, b_end);
    lp.b = Eigen::VectorXd::Map(bv.data(), static_cast<Eigen::Index>(bv.size()));

    // Parse c
    std::size_t c_start = find_value_start(json, "c");
    std::size_t c_end;
    auto cv = parse_number_array(json, c_start, c_end);
    lp.c = Eigen::VectorXd::Map(cv.data(), static_cast<Eigen::Index>(cv.size()));

    return lp;
}

// ---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: bench_mps <data.json> [N]\n";
        return 1;
    }
    std::string json_path = argv[1];
    int n_runs = (argc >= 3) ? std::atoi(argv[2]) : 10;

    LpData data = load_lp(json_path);
    LinearProgram lp_proto(data.A, data.b, data.c);

    using Clock = std::chrono::high_resolution_clock;
    using Dur   = std::chrono::duration<double, std::milli>;

    double last_z = 0.0;
    int    last_iterations = 0;
    std::vector<double> times_ms;
    times_ms.reserve(n_runs);

    for (int r = 0; r < n_runs; ++r) {
        LinearProgram lp(lp_proto);  // copy so each solve gets fresh data
        auto t0 = Clock::now();
        SimplexResult res = solve(lp);
        auto t1 = Clock::now();
        times_ms.push_back(Dur(t1 - t0).count());
        last_z          = res.z;
        last_iterations = res.iterations;
    }

    // Emit JSON
    std::cout << "{\"z\": " << std::fixed;
    std::cout.precision(4);
    std::cout << last_z;
    std::cout << ", \"iterations\": " << last_iterations;
    std::cout << ", \"times_ms\": [";
    std::cout.precision(6);
    for (int i = 0; i < n_runs; ++i) {
        if (i) std::cout << ", ";
        std::cout << times_ms[i];
    }
    std::cout << "]}\n";
    return 0;
}
