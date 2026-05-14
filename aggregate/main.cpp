#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Table
{
    std::vector<std::string> header;
    std::vector<std::vector<long double>> rows;
};

static FILE* open_file_to_write(const std::filesystem::path& filenamePath)
{
    std::string filename = filenamePath.string();
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open write file " << filename << std::endl;
        exit(1);
    }
    return file;
}

static std::string trim(const std::string& s)
{
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

static bool is_comment_line(const std::string& line)
{
    const std::string t = trim(line);
    return t.empty() || t[0] == '#';
}

static Table read_table(const std::filesystem::path& path)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open read file " + path.string());
    }

    Table out;
    std::string line;
    while (std::getline(file, line)) {
        if (is_comment_line(line)) {
            out.header.push_back(line);
            continue;
        }

        std::stringstream ss(line);
        std::vector<long double> row;
        long double value;
        while (ss >> value) {
            row.push_back(value);
        }
        if (!row.empty()) {
            out.rows.push_back(std::move(row));
        }
    }

    return out;
}

static bool parse_range_header(const std::vector<std::string>& header, long int& start, long int& stop)
{
    if (header.empty()) return false;

    std::stringstream ss(trim(header.front()));
    char hash = '\0';
    if (!(ss >> hash)) return false;
    if (hash != '#') return false;
    if (!(ss >> start >> stop)) return false;
    return true;
}

static long double mean(const std::vector<long double>& v)
{
    if (v.empty()) return 0.0L;
    long double s = 0.0L;
    for (long double x : v) s += x;
    return s / static_cast<long double>(v.size());
}

static long double second_moment(const std::vector<long double>& v)
{
    if (v.empty()) return 0.0L;
    long double s = 0.0L;
    for (long double x : v) s += x * x;
    return s / static_cast<long double>(v.size());
}

static std::filesystem::path make_input_path(const std::string& stem, int L, int T, int num)
{
    std::ostringstream ss;
    ss << "../res/" << stem << "_L" << L << "_T" << T << "_num" << num << ".txt";
    return ss.str();
}

static std::filesystem::path make_output_path(const std::string& stem, int L, int T, int first_num, int last_num, const std::string& mode, int col = -1)
{
    std::ostringstream ss;
    ss << "../agg_res/" << stem << "_L" << L << "_T" << T << "_num" << first_num << "_to_" << last_num;
    if (col >= 0) ss << "_col" << col;
    ss << "_" << mode << ".txt";
    return ss.str();
}

static void ensure_output_dir()
{
    std::filesystem::create_directories("../agg_res");
}

static void aggregate_curve(const std::string& stem, int L, int T, int first_num, int last_num)
{
    const int nfiles = last_num - first_num + 1;
    if (nfiles <= 0) throw std::runtime_error("Invalid run range");

    std::vector<std::vector<long double>> sums;
    std::size_t nrows = 0;
    std::size_t ncols = 0;
    long int range_start = 0;
    long int range_stop = 0;
    bool have_range = false;
    int found_files = 0;

    for (int num = first_num; num <= last_num; ++num) {
        const auto path = make_input_path(stem, L, T, num);
        if (!std::filesystem::exists(path)) continue;
        const Table table = read_table(path);
        if (table.rows.empty()) {
            throw std::runtime_error("No numeric rows found in " + path.string());
        }

        if (!have_range) {
            have_range = parse_range_header(table.header, range_start, range_stop);
        }
        if (sums.empty()) {
            nrows = table.rows.size();
            ncols = table.rows.front().size();
            sums.assign(nrows, std::vector<long double>(ncols, 0.0L));
        } else {
            if (table.rows.size() != nrows) {
                throw std::runtime_error("Mismatching row count in " + path.string());
            }
            if (table.rows.front().size() != ncols) {
                throw std::runtime_error("Mismatching column count in " + path.string());
            }
        }

        for (std::size_t i = 0; i < nrows; ++i) {
            if (table.rows[i].size() != ncols) {
                throw std::runtime_error("Mismatching column count on row " + std::to_string(i) + " in " + path.string());
            }
            for (std::size_t j = 0; j < ncols; ++j) {
                sums[i][j] += table.rows[i][j];
            }
        }
        ++found_files;
    }

    if (found_files == 0) throw std::runtime_error("No input files found in requested range");

    const auto out_path = make_output_path(stem, L, T, first_num, last_num, "avg");
    FILE* out = open_file_to_write(out_path);
    if (have_range) {
        fprintf(out, "# %ld %ld\n", range_start, range_stop);
    }
    fprintf(out, "# Aggregated over runs\n");

    for (std::size_t i = 0; i < nrows; ++i) {
        for (std::size_t j = 0; j < ncols; ++j) {
            sums[i][j] /= static_cast<long double>(found_files);
        }
        for (std::size_t j = 0; j < ncols; ++j) {
            if (j + 1 == ncols) fprintf(out, "%.12Lf", sums[i][j]);
            else fprintf(out, "%.12Lf ", sums[i][j]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
}

static std::vector<long double> collect_samples(const std::string& stem, int L, int T, int first_num, int last_num, int col)
{
    if (col < 0) throw std::runtime_error("Column index required for sample-based aggregation");

    std::vector<long double> values;
    for (int num = first_num; num <= last_num; ++num) {
        const auto path = make_input_path(stem, L, T, num);
        if (!std::filesystem::exists(path)) continue;
        const Table table = read_table(path);
        for (const auto& row : table.rows) {
            if (static_cast<std::size_t>(col) >= row.size()) {
                throw std::runtime_error("Requested column out of range in " + path.string());
            }
            values.push_back(row[col]);
        }
    }
    return values;
}

static void aggregate_samples(const std::string& stem, int L, int T, int first_num, int last_num, int col)
{
    std::vector<long double> values = collect_samples(stem, L, T, first_num, last_num, col);
    if (values.empty()) {
        throw std::runtime_error("No values collected for sample aggregation");
    }

    const long double m1 = mean(values);
    const long double m2 = second_moment(values);
    const auto out_path = make_output_path(stem, L, T, first_num, last_num, "samples", col);
    FILE* out = open_file_to_write(out_path);
    fprintf(out, "# %.12Lf %.12Lf\n", m1, m2);
    for (long double v : values) {
        fprintf(out, "%.12Lf\n", v);
    }
    fclose(out);
}

static void aggregate_hist(const std::string& stem, int L, int T, int first_num, int last_num, int col)
{
    std::vector<long double> values = collect_samples(stem, L, T, first_num, last_num, col);
    if (values.empty()) {
        throw std::runtime_error("No values collected for histogram aggregation");
    }

    std::sort(values.begin(), values.end());
    const long double m1 = mean(values);
    const long double m2 = second_moment(values);

    std::vector<long double> bins;
    std::vector<unsigned long long> counts;

    for (long double v : values) {
        if (bins.empty() || bins.back() != v) {
            bins.push_back(v);
            counts.push_back(1);
        } else {
            ++counts.back();
        }
    }

    const auto out_path = make_output_path(stem, L, T, first_num, last_num, "hist", col);
    FILE* out = open_file_to_write(out_path);
    fprintf(out, "# %.12Lf %.12Lf\n", m1, m2);
    for (std::size_t i = 0; i < bins.size(); ++i) {
        fprintf(out, "%.12Lf %llu\n", bins[i], counts[i]);
    }
    fclose(out);
}

static void usage(const char* name)
{
    std::cerr << "Usage: " << name
              << " <L> <T> <FIRST_NUM> <LAST_NUM> <STEM> <MODE> [COLUMN]\n"
              << "Modes: curve, samples, hist\n";
}

int main(int argc, char* argv[])
{
    if (argc < 7) {
        usage(argv[0]);
        return 1;
    }

    const int L = std::atoi(argv[1]);
    const int T = std::atoi(argv[2]);
    const int first_num = std::atoi(argv[3]);
    const int last_num = std::atoi(argv[4]);
    const std::string stem = argv[5];
    const std::string mode = argv[6];
    const int col = (argc >= 8) ? std::atoi(argv[7]) : -1;

    ensure_output_dir();

    if (mode == "curve") {
        aggregate_curve(stem, L, T, first_num, last_num);
    } else if (mode == "samples") {
        aggregate_samples(stem, L, T, first_num, last_num, col);
    } else if (mode == "hist") {
        aggregate_hist(stem, L, T, first_num, last_num, col);
    } else {
        throw std::runtime_error("Unknown mode: " + mode);
    }

    return 0;
}
