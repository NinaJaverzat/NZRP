/*
Full-range convolutor.

This is the generic full-range convolution tool. The binomial logic is copied
from the existing implementation; only the I/O contract is generalized so it can
read different raw file families.
*/

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <stack>
#include <vector>

#include <math.h>
#include <chrono>

#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <sys/resource.h>

#include <bits/stdc++.h>
#include <boost/math/distributions/binomial.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <initializer_list>

struct Table
{
    std::vector<std::string> header;
    std::vector<std::vector<long double>> rows;
};

FILE* open_file_to_write(const std::ostringstream& filenameStream)
{
    std::string filename = filenameStream.str();
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

static Table read_table(const std::ostringstream& filenameStream)
{
    std::string filename = filenameStream.str();
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open read file " + filename);

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
        if (!row.empty()) out.rows.push_back(std::move(row));
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

static std::vector<long double> read_column_curve(const std::ostringstream& filenameStream, const int& col, const std::size_t& num_rows)
{
    std::string filename = filenameStream.str();
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open read file " + filename);

    std::vector<long double> vec(num_rows, 0.0L);
    std::string line;
    std::size_t i = 0;
    while (std::getline(file, line)) {
        if (is_comment_line(line)) continue;
        std::stringstream ss(line);
        std::vector<long double> row;
        long double value;
        while (ss >> value) row.push_back(value);
        if (row.empty()) continue;
        if (static_cast<std::size_t>(col) >= row.size()) {
            throw std::runtime_error("Requested column out of range in " + filename);
        }
        if (i >= num_rows) throw std::runtime_error("Unexpected number of rows in " + filename);
        vec[i] = row[col];
        ++i;
    }
    if (i != num_rows) throw std::runtime_error("Unexpected number of rows in " + filename);
    return vec;
}

static long double get_RL_val(const std::vector<size_t>& m_vals, const std::vector<long double>& cdf,
                              const size_t& m_min, const size_t& m_max, const size_t m_val)
{
    long double RL_val;

    if (m_val < m_min) RL_val = 0.0;
    else if (m_val > m_max) RL_val = 1.0;
    else {
        auto it = std::lower_bound(m_vals.begin(), m_vals.end(), m_val);
        size_t RL_ind = it - m_vals.begin();
        if (*it == m_val) RL_val = cdf[RL_ind];
        if (*it > m_val) RL_val = cdf[RL_ind-1];
        if (*it < m_val) std::cout << "UNEXPECTED INDEX FOUND\n";
    }

    return RL_val;
}

static std::vector<size_t> cdf_from_samples(std::vector<size_t> data, std::vector<long double>& cdf_out)
{
    std::sort(data.begin(), data.end());
    std::vector<size_t> vals;
    if (data.empty()) return vals;

    vals.push_back(static_cast<size_t>(data[0]));
    cdf_out.push_back(1.0L);
    for (std::size_t i = 1; i < data.size(); ++i) {
        const size_t m = data[i];
        if (m == vals.back()) ++cdf_out.back();
        else {
            vals.push_back(m);
            cdf_out.push_back(cdf_out.back() + 1.0L);
        }
    }
    const long double norm = static_cast<long double>(data.size());
    for (auto& x : cdf_out) x /= norm;
    return vals;
}

int main(int argc, char* argv[])
{
    if (argc < 8) {
        std::cerr << "Usage: " << argv[0]
                  << " <L> <T> <FIRST_NUM> <LAST_NUM> <STEM> <MODE> <NCOLS> [COLS...] [MODEL]\n"
                  << "Modes: pair, curve, cdf, wrap\n";
        return 1;
    }

    const int L = atoi(argv[1]);
    const int T = atoi(argv[2]);
    const int first_num = atoi(argv[3]);
    const int last_num = atoi(argv[4]);
    const std::string stem = argv[5];
    const std::string mode = argv[6];
    const int ncols = atoi(argv[7]);
    if (ncols <= 0) throw std::runtime_error("NCOLS must be positive");
    if (argc != 8 + ncols && argc != 9 + ncols) throw std::runtime_error("Wrong number of columns passed");

    std::vector<int> cols(ncols, 0);
    for (int i = 0; i < ncols; ++i) cols[i] = atoi(argv[8 + i]);
    const std::string model = (argc == 9 + ncols) ? argv[8 + ncols] : stem;

    const int N = L*L;
    const int M = 3*N;
    const int nfiles = last_num - first_num + 1;
    if (nfiles <= 0) throw std::runtime_error("Invalid run range");
    const std::string outstem = (mode == "cdf" || mode == "wrap") ? model : stem;

    std::ostringstream fname;

    std::vector<std::vector<long double>> curves(ncols, std::vector<long double>(M, 0.0L));
    std::vector<std::vector<long double>> convs(ncols, std::vector<long double>(M, 0.0L));
    std::vector<long double> binomial(M + 1, 0.0L);
    int found_files = 0;

    fname << "../conv_res/" << outstem << "_L" << L << "_T" << T << "_num" << first_num << "_to_" << last_num << "_convoluted.txt";
    FILE* myOUTfile = open_file_to_write(fname);
    fname.str("");
    fname.clear();

    if (mode == "pair" || mode == "curve") {
        for (int i = first_num; i <= last_num; ++i) {
            fname << "../res/" << stem << "_L" << L << "_T" << T << "_num" << i << ".txt";
            if (!std::filesystem::exists(fname.str())) {
                fname.str("");
                fname.clear();
                continue;
            }
            auto table = read_table(fname);
            if (table.rows.size() != static_cast<std::size_t>(M)) {
                throw std::runtime_error("Unexpected row count in " + fname.str());
            }
            for (int c = 0; c < ncols; ++c) {
                for (int m = 0; m < M; ++m) {
                    if (static_cast<std::size_t>(cols[c]) >= table.rows[m].size()) {
                        throw std::runtime_error("Column out of range in " + fname.str());
                    }
                    curves[c][m] += table.rows[m][cols[c]];
                }
            }
            fname.str("");
            fname.clear();
            ++found_files;
        }

        if (found_files == 0) throw std::runtime_error("No input files found in requested range");

        const long double norm = static_cast<long double>(T) * static_cast<long double>(found_files);
        for (int c = 0; c < ncols; ++c) {
            std::transform(curves[c].begin(), curves[c].end(), curves[c].begin(), [norm](auto& x) { return x / norm; });
        }

        long double value;
        long double dp = 1.0 / M;
        long double p = dp;
        int n;
        int first_index = 0;
        int last_index = 0;

        for (int i = 0; i < M - 1; ++i) {
            boost::math::binomial_distribution<long double> dist(M, p);
            for (n = first_index + 1; n <= M; ++n) {
                const long double val = boost::math::pdf(dist, n);
                if (val > 1e-12) break;
            }
            first_index = n - 1;
            for (n = first_index + 1; n <= M; ++n) {
                const long double val = boost::math::pdf(dist, n);
                binomial[n - 1] = val;
                if (val < 1e-12 || n == M) break;
            }
            last_index = n;

            fprintf(myOUTfile, "%.12Lf", p);
            for (int c = 0; c < ncols; ++c) {
                value = 0;
                for (n = first_index; n <= last_index; ++n) value += curves[c][n] * binomial[n];
                convs[c][i] = value;
                fprintf(myOUTfile, " %.12Lf %.12Lf", curves[c][i], value);
            }
            fprintf(myOUTfile, "\n");
            p += dp;
        }

        boost::math::binomial_distribution<long double> dist(M, 1.0);
        for (n = first_index + 1; n <= M; ++n) {
            const long double val = boost::math::pdf(dist, n);
            if (val > 1e-12) break;
        }
        first_index = n - 1;
        for (n = first_index + 1; n <= M; ++n) {
            const long double val = boost::math::pdf(dist, n);
            binomial[n - 1] = val;
            if (val < 1e-12 || n == M) break;
        }
        last_index = n;

        fprintf(myOUTfile, "%.12Lf", 1.0L);
        for (int c = 0; c < ncols; ++c) {
            value = 0;
            for (n = first_index; n <= last_index; ++n) value += curves[c][n] * binomial[n];
            convs[c][M - 1] = value;
            fprintf(myOUTfile, " %.12Lf %.12Lf", curves[c][M - 1], value);
        }
        fprintf(myOUTfile, "\n");
    } else if (mode == "cdf" || mode == "wrap") {
        if (ncols != 2) throw std::runtime_error("cdf mode expects exactly two columns");
        std::vector<size_t> xraw, xyraw;
        for (int i = first_num; i <= last_num; ++i) {
            fname << "../res/" << stem << "_L" << L << "_T" << T << "_num" << i << ".txt";
            if (!std::filesystem::exists(fname.str())) {
                fname.str("");
                fname.clear();
                continue;
            }
            auto table = read_table(fname);
            for (const auto& row : table.rows) {
                if (static_cast<std::size_t>(cols[0]) >= row.size() || static_cast<std::size_t>(cols[1]) >= row.size()) {
                    throw std::runtime_error("Requested column out of range in " + fname.str());
                }
                xraw.push_back(static_cast<size_t>(row[cols[0]]));
                xyraw.push_back(static_cast<size_t>(row[cols[1]]));
            }
            fname.str("");
            fname.clear();
            ++found_files;
        }

        if (found_files == 0) throw std::runtime_error("No input files found in requested range");

        std::vector<long double> xcdf, xycdf;
        auto xvals = cdf_from_samples(xraw, xcdf);
        auto xyvals = cdf_from_samples(xyraw, xycdf);

        if (xvals.empty() || xyvals.empty()) throw std::runtime_error("Empty data in cdf mode");

        size_t end = (xcdf.size() < xycdf.size()) ? xcdf.size() : xycdf.size();

        std::ostringstream avgname;
        avgname << "../conv_res/" << model << "R_L" << L << "_T" << T << "_num" << first_num << "_to_" << last_num << "_averaged.txt";
        FILE* myAVGfile = open_file_to_write(avgname);

        fprintf(myAVGfile, "# Columns: x_threshold x_cdf xy_threshold xy_cdf\n");
        fprintf(myAVGfile, "%lu %.16Lf %lu %.16Lf\n", xvals[0] - 1, 0.0L, xyvals[0] - 1, 0.0L);
        for (size_t i = 0; i < end; ++i) {
            fprintf(myAVGfile, "%lu %.16Lf %lu %.16Lf\n", xvals[i], xcdf[i], xyvals[i], xycdf[i]);
        }
        if (end == xcdf.size()) {
            for (size_t i = end; i < xycdf.size(); ++i) {
                fprintf(myAVGfile, "%lu %.16Lf %lu %.16Lf\n", xvals[end-1] + (i + 1 - end), 1.0L, xyvals[i], xycdf[i]);
            }
        } else if (end == xycdf.size()) {
            for (size_t i = end; i < xcdf.size(); ++i) {
                fprintf(myAVGfile, "%lu %.16Lf %lu %.16Lf\n", xvals[i], xcdf[i], xyvals[end-1] + (i + 1 - end), 1.0L);
            }
        }
        fclose(myAVGfile);

        size_t m1 = xvals[0];
        size_t m2 = xyvals.back();
        if (xyvals[0] < m1) {
            m1 = xyvals[0];
            std::cout << "The smallest mxy results smaller than the smallest mx. This is unexpected!" << std::endl;
        }
        if (xvals.back() > m2) {
            m2 = xvals.back();
            std::cout << "The largest mx results larger than the largest mxy. This is unexpected!" << std::endl;
        }

        long double dp = 1.0L / static_cast<long double>(M);
        long double p = static_cast<long double>(m1) / static_cast<long double>(M);

        while (true) {
            boost::math::binomial_distribution<long double> dist(M, p);
            const long double val = boost::math::pdf(dist, static_cast<long int>(m1));
            if (val < 1e-12L) break;
            p -= dp;
        }

        long double p0 = p;
        size_t last_m = m1;
        size_t first_m = m1;

        bool mass_found = 0;
        while (true) {
            first_m -= 1;
            boost::math::binomial_distribution<long double> dist(M, p0);
            const long double val = boost::math::pdf(dist, static_cast<long int>(first_m));
            if (mass_found == 0 && val > 1e-11L) mass_found = 1;
            if (mass_found && val < 1e-12L) break;
        }

        p0 += dp;
        for (p = p0; p <= 1.0L; p += dp) {
            while (true) {
                boost::math::binomial_distribution<long double> dist(M, p);
                const long double val = boost::math::pdf(dist, static_cast<long int>(first_m));
                if (val > 1e-12L || first_m > m2) break;
                ++first_m;
            }

            last_m = first_m;
            std::vector<long double> binomial_local = {};
            while (true) {
                boost::math::binomial_distribution<long double> dist(M, p);
                const long double val = boost::math::pdf(dist, static_cast<long int>(last_m));
                binomial_local.push_back(val);
                if (val < 1e-12L || last_m == static_cast<size_t>(M)) break;
                ++last_m;
            }

            std::reverse(binomial_local.begin(), binomial_local.end());

            long double RL_x_conv = 0.0L;
            long double RL_xy_conv = 0.0L;
            for (size_t bin_ind = 0; bin_ind < binomial_local.size(); ++bin_ind) {
                const long double bin_val = binomial_local[bin_ind];
                const size_t m_val = first_m + bin_ind;

                const long double RL_x = get_RL_val(xvals, xcdf, xvals[0], xvals.back(), m_val);
                const long double RL_xy = get_RL_val(xyvals, xycdf, xyvals[0], xyvals.back(), m_val);

                RL_x_conv += bin_val * RL_x;
                RL_xy_conv += bin_val * RL_xy;
            }

            fprintf(myOUTfile, "%.16Lf %.16Lf %.16Lf %.16Lf %.16Lf\n",
                    p, RL_x_conv, RL_xy_conv, 2.0L * RL_x_conv - RL_xy_conv, RL_x_conv - RL_xy_conv);
            fflush(myOUTfile);

            if (first_m > m2) break;
        }
    } else {
        throw std::runtime_error("Unknown mode: " + mode);
    }

    fclose(myOUTfile);
    return 0;
}
