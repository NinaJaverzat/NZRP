/*
Convolutor.

This is the partial-range convolution tool. The binomial logic is copied
from the existing implementation; only the I/O contract is generalized so it can
read different raw file families and explicit run ranges.
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

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

static FILE* open_file_to_write(const std::ostringstream& filenameStream)
{
    std::string filename = filenameStream.str();
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open write file " << filename << std::endl;
        exit(1);
    }
    return file;
}

static std::pair<size_t, size_t> read_range(const std::ostringstream& filenameStream)
{
    std::string filename = filenameStream.str();
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open read file " + filename);

    std::string hash;
    size_t m1, m2;
    file >> hash >> m1 >> m2;

    if (hash != "#") throw std::runtime_error("Invalid header format in " + filename);
    return {m1, m2};
}

static std::string first_existing_input(const std::string& model, const int L, const int T, const int first_num, const int last_num)
{
    std::ostringstream fname;
    for (int i = first_num; i <= last_num; ++i)
    {
        fname << "../res/" << model << "_L" << L << "_T" << T << "_num" << i << ".txt";
        if (std::filesystem::exists(fname.str())) return fname.str();
        fname.str("");
        fname.clear();
    }
    throw std::runtime_error("No input files found in requested range");
}

static void read_input(const std::ostringstream& filenameStream,
                       std::vector<long double>& Pinf_m,
                       std::vector<long int>& Pinf_norm,
                       std::vector<long double>& Chi_m,
                       const int& num_lines,
                       const int pinf_col,
                       const int norm_col,
                       const int chi_col)
{
    std::string filename = filenameStream.str();
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Failed to open read file " + filename);

    std::string line;
    bool header_skipped = false;
    int i = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string first;
        ss >> first;
        if (first.empty()) continue;
        if (first[0] == '#') {
            continue;
        }

        std::vector<long double> row;
        long double value;
        std::stringstream rowss(line);
        while (rowss >> value) row.push_back(value);
        if (row.empty()) continue;

        const std::size_t pcol = static_cast<std::size_t>(pinf_col);
        const std::size_t ncol = static_cast<std::size_t>(norm_col);
        const std::size_t ccol = static_cast<std::size_t>(chi_col);
        if (pcol >= row.size() || ncol >= row.size() || ccol >= row.size()) {
            throw std::runtime_error("Requested column out of range in " + filename);
        }

        if (!header_skipped) {
            header_skipped = true;
        }

        Pinf_m[i] += row[pcol];
        Pinf_norm[i] += static_cast<long int>(row[ncol]);
        Chi_m[i] += row[ccol];
        ++i;
    }

    if (i != num_lines) {
        throw std::runtime_error("Unexpected number of lines: " + std::to_string(i) + " instead of " + std::to_string(num_lines));
    }
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// MAIN FUNCTION

int main(int argc, char* argv[]) {

    if (argc != 7 && argc != 9) {
      std::cerr << "Usage: " << argv[0]
                << " <LATTICE LINEAR SIZE> <NUMBER OF SAMPLES> <FIRST RUN> <LAST RUN> <STEM> <PINF COL> [NORM COL CHI COL]"
                << std::endl;
      return 1;
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    const int L = atoi(argv[1]);
    const int T = atoi(argv[2]);
    const int first_num = atoi(argv[3]);
    const int last_num = atoi(argv[4]);
    const std::string model = argv[5];
    const int pinf_col = atoi(argv[6]);
    const int norm_col = (argc == 9) ? atoi(argv[7]) : 1;
    const int chi_col = (argc == 9) ? atoi(argv[8]) : 2;

    const int N = L*L;
    const int M = 3*N;
    const int num_files = last_num - first_num + 1;
    if (num_files <= 0) throw std::runtime_error("Invalid run range");
    int found_files = 0;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// OUTPUT FILES TO BE WRITTEN

    std::ostringstream fname;

    fname << "../conv_res/" << model << "_L" << L << "_T" << T << "_num" << first_num << "_to_" << last_num << "_averaged.txt";
    FILE* myOUTfile_avg = open_file_to_write(fname);

    fname.str("");
    fname.clear();

    fname << "../conv_res/" << model << "_L" << L << "_T" << T << "_num" << first_num << "_to_" << last_num << "_convoluted.txt";
    FILE* myOUTfile_conv = open_file_to_write(fname);

    fname.str("");
    fname.clear();

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// READ THE RANGE AND THEN THE INPUT

    fname << first_existing_input(model, L, T, first_num, last_num);
    auto [m1, m2] = read_range(fname);
    fname.str("");
    fname.clear();

    fprintf(myOUTfile_avg, "# %lu %lu\n", m1, m2);

    const size_t len = m2 - m1 + 1;
    std::vector<long int> Pinf_norm(len, 0);
    std::vector<long double> Pinf_m(len, 0), Chi_m(len, 0);
    std::vector<long double> Pinf_conv, Chi_conv;

    for (int i = first_num; i <= last_num; ++i)
    {
        fname << "../res/" << model << "_L" << L << "_T" << T << "_num" << i << ".txt";
        if (!std::filesystem::exists(fname.str())) {
            fname.str("");
            fname.clear();
            continue;
        }
        read_input(fname, Pinf_m, Pinf_norm, Chi_m, static_cast<int>(len), pinf_col, norm_col, chi_col);
        fname.str("");
        fname.clear();
        ++found_files;
    }

    if (found_files == 0) throw std::runtime_error("No input files found in requested range");

// NORMALIZE THE INPUT

    std::transform(Pinf_m.begin(), Pinf_m.end(), Pinf_norm.begin(), Pinf_m.begin(), [](long double a, long double b) { return a / b; });
    Pinf_norm.clear();
    Pinf_norm.shrink_to_fit();

    const long double norm = static_cast<long double>(found_files) * static_cast<long double>(T == -1 ? 1 : T);
    std::transform(Chi_m.begin(), Chi_m.end(), Chi_m.begin(), [norm](auto& x) { return x / norm; });

// CONVOLUTE

    /*
    Part 1: find first relevant value of p.
    The smallest value of p is such that the binomial is zero for the last time every m <= m1.
    Given that P_inf = 0 at m = m1, we know that sum_m B(m, M, p) * P_inf(m) = 0 at the chosen p:
    - Where B(m, M, p) > 1e-12, P_inf = 0
    - Where B(m, M, p) < 1e-12, the summation is taken as zero by construction
    */

    long double pstar;
    long double dp = 1.0 / M;
    long double p = static_cast<long double>(m1) / static_cast<long double>(M);

    while (true)
    {
        boost::math::binomial_distribution<long double> dist(M, p);
        const long double val = boost::math::pdf(dist, m1);
        if (val < 1e-12L) break;
        p -= dp;
    }

    long double p0 = p;

    /*
    Part 2: find the domain of the first binomial.
    At this stage, we know that at p = p0 the binomial has basically zero weigth for m = m1 so that the convolution is zero at this p.
    Here we find the domain of the binomail where it is non-zero at the given p. All the summations at future values of p will be made
    by first computing the new domain given the new p and then computing the sum. The new domain will be computed by shifting the previous one.
    */

    size_t last_m = m1;
    size_t first_m = m1;

    bool mass_found = 0;
    while (true)
    {
        first_m -= 1;
        boost::math::binomial_distribution<long double> dist(M, p0);
        const long double val = boost::math::pdf(dist, first_m);
        if (mass_found == 0 && val > 1e-11L) mass_found = 1;
        if (mass_found && val < 1e-12L) break;
    }

    /*
    At this point we have built a binomial distribution with parameter p0, B(m, M, p0) and we know that:
    - B(m, M, p0) < 1e-12 ~ 0 for m < first_m
    - B(m, M, p0) < 1e-12 ~ 0 for m > last_m
    - last_m < m1
    - P_inf (m = m1) = 0

    So that

    sum_m B(m, M, p0) * P_inf(m) = 0

    We can now perform

    PART 3: the convolution.

    We proceed increasing p until we reach a value of p that is so large that the
    binomial is non-zero if m > m2, which is the last value for which we have measured P_inf.
    In this way we convolute the percolation strength up to the value of p that our measure allows to assess correctly
    The first value of p that we can not assess correctly is a p such that we have

    B(m2+1, M, p) > 1e-12 ~ 0

    Indeed, at that value of p we would need to know P_inf(m2+1), which we do not know. So we can not compute P_inf(p).
    We quit when such p is reached.

    This is the convolution.
    */

    pstar = -1;
    long double Pinf_val, Chi_val;

    p0 += dp;
    for (p = p0; p <= 1; p += dp)
    {
        boost::math::binomial_distribution<long double> dist(M, p);
        while (true)
        {
            const long double val = boost::math::pdf(dist, first_m);
            if (val > 1e-12L) break;
            ++first_m;
        }

        last_m = first_m;
        std::vector<long double> binomial_local = {};
        while (true)
        {
            const long double val = boost::math::pdf(dist, last_m);
            binomial_local.push_back(val);
            if (val < 1e-12L || last_m > m2) break;
            ++last_m;
        }

        if (last_m > m2) break;

        std::reverse(binomial_local.begin(), binomial_local.end());

        Pinf_val = Chi_val = 0;
        for (size_t bin_ind = 0; bin_ind < binomial_local.size(); ++bin_ind)
        {
            const long double bin_val = binomial_local[bin_ind];
            size_t arr_ind = (first_m + bin_ind > m1) ? first_m + bin_ind - m1 : 0;

            if (arr_ind > len)
            {
                std::cout << "ERROR! arr_ind = " << arr_ind << " but len = " << len << std::endl;
                std::cout << "Note that (m1, m2) = (" << m1 << ", " << m2 << ")" << std::endl;
                std::cout << "And that (first_m, last_m) = (" << first_m << ", " << last_m << ")" << std::endl;
                std::cout << "The error happens at bin_ind = " << bin_ind << std::endl;
                exit(1);
            }

            Pinf_val += bin_val * Pinf_m[arr_ind];
            Chi_val += bin_val * Chi_m[arr_ind];
        }
        Pinf_conv.push_back(Pinf_val);
        Chi_conv.push_back(Chi_val);

        if (first_m >= m1 && pstar < 0) pstar = p;
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// PRINT OUTPUT

    for (size_t i = 0; i < len; ++i) {
        fprintf(myOUTfile_avg, "%.16Lf %.16Lf\n", Pinf_m[i], Chi_m[i]);
        fflush(myOUTfile_avg);
    }

    fprintf(myOUTfile_conv, "# %.16Lf %.16Lf %.16Lf\n", p0, p, pstar);
    for (size_t i = 0; i < Pinf_conv.size(); ++i) {
        fprintf(myOUTfile_conv, "%.16Lf %.16Lf\n", Pinf_conv[i], Chi_conv[i]);
        fflush(myOUTfile_conv);
    }

    return 0;
}
