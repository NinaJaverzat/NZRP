/*
Convolutor.

This code takes in input L, T and a third variable, just like NZ.
The third variable represents the number of parallel simulations made of the NZ code using L and T
So, if we make NZ simulations from num=0 to num=10, this third number will be 10. Let's call this number TOT.

Given these information, this code reads TOT result files of the NZ code, which are always written as a function of the bond number m
and makes the convolution with the binomial distribution, producing a single file contains
1) The average over TOT files of the results as a function of m
2) The convoluted average as a function of p

This is done for all the quantities of interest (so COP and its fluctuations, ROP and its fluctuations, Time per bond and, if we want
other observables like the number of searches and the number of ended searches)
*/

#include <cstdio>
#include <fstream>
#include <iostream>

#include <stack>
#include <vector>

#include <math.h>
// #include <time.h>				// Get clock time, to measure total run time
#include <chrono>           // Better clock time

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

void read_input(const std::ostringstream& filenameStream, std::vector<long double>& COP, std::vector<long double>& ROP, const int& num_lines)
{

  std::string filename = filenameStream.str();
  std::ifstream file(filename);
  if (!file.is_open()) if (!file.is_open()) throw std::runtime_error("Failed to open read file " + filename);

  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Skip first two lines
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  int i = 0;
  long double vals[2];
  while (file >> vals[0] >> vals[1])
  {
      COP[i] += vals[0];
      ROP[i] += vals[1];
      ++i;
  }

  if (i != num_lines) throw std::runtime_error("Unexpected number of lines: " + std::to_string(i) +  " instead of " + std::to_string(num_lines));

}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


// MAIN FUNCTION

int main(int argc, char* argv[]) {

    if (argc!=4) {
      std::cerr << "Usage: " << argv[0] << " <LATTICE LINEAR SIZE>  <NUMBER OF SAMPLES> <NUMBER OF RUNS> " << std::endl;
      return 1;
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


    const int L = atoi(argv[1]);                                                // Linear size of the lattice
    const int T = atoi(argv[2]);                                                // Number of trials
    const int TOT = atoi(argv[3]);                                              // Number of runs

    const int N = L*L;
    const int M = 3*N;

    std::vector<long double> COP(M, 0), ROP(M, 0);
    std::vector<long double> COP_conv(M, 0), ROP_conv(M, 0);
    std::vector<long double> binomial(M, 0);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// OUTPUT FILES TO BE WRITTEN

    std::ostringstream fname;

    fname << "../conv_res/Smax_L" << L << "_T" << T << "_TOT" << TOT << "_convoluted.txt";               // Connectivity order parameter vs number of active bonds
    FILE* myOUTfile = open_file_to_write(fname);
    fprintf(myOUTfile, "# Columns: averge Smax in CP, convluted average Smax in CP, averge Smax in RP, convluted average Smax in RP\n");

    fname.str("");
    fname.clear();

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// READ THE INPUT

    for (int i=1; i<=TOT; ++i)                                                              // NOTE THAT THE NUMBER OF FILES STATS AT 1 AND ENDS AT TOT
    {
        fname << "../res/Smax_L" << L << "_T" << T << "_num" << i << ".txt";

        read_input(fname, COP, ROP, M);

        fname.str("");
        fname.clear();
    }

    // print_vecs(COP1, COP2);

// NORMALIZE THE INPUT

    double norm = T * (TOT+1);
    std::transform(COP.begin(), COP.end(), COP.begin(), [norm](auto& x) { return x / norm; });
    std::transform(ROP.begin(), ROP.end(), ROP.begin(), [norm](auto& x) { return x / norm; });

// CONVOLUTE


    long double value;
    long double dp = 1.0 / M;
    long double p = dp;


    int n;
    int first_index = 0;
    int last_index = 0;

    for (int i=0; i<M-1; ++i)
    {
        boost::math::binomial_distribution<long double> dist(M, p);


        // Find the first value at which the binomial is larger than 1e-12 given the new value of p
        // We start from p = 1/M and n = 1 successes
        // Note that the value of n is the first relevant index + 1
        for (n=first_index+1; n<=M; ++n)
        {
            const long double val = boost::math::pdf(dist, n);
            if (val > 1e-12) break;
        }

        first_index = n-1;                                                      // We start convoluting at this index, i.e., we give for granted that previous addednds where 0 as the binomial is too small (given the precision of long double)

        // Find the last value at which the binomial is larger than 1e-12 given the new value of p
        // We start from p = 1/M and n = 1 successes
        // Note that the value of n is also the last relevant index
        for (n=first_index+1; n<=M; ++n)
        {
            const long double val = boost::math::pdf(dist, n);
            binomial[n-1] = val;
            if (val < 1e-12 || n == M) break;
        }
        last_index = n;

        // std::cout << "For i = " << i << " I have p = " << p << " I find first_index = " << first_index << " and last_index = " << last_index << std::endl;

        value = 0;
        for (n=first_index; n<=last_index; n++) { value += COP[n]*binomial[n]; }
        COP_conv[i] = value;

        value = 0;
        for (n=first_index; n<=last_index; n++) { value += ROP[n]*binomial[n]; }
        ROP_conv[i] = value;

        p += dp;

    }


    // We do the last iteration manually as p might become > 1 due to floating point problems

    boost::math::binomial_distribution<long double> dist(M, 1.0);


    // Find the first value at which the binomial is larger than 1e-12 given the new value of p
    // We start from p = 1/M and n = 1 successes
    // Note that the value of n is the first relevant index + 1
    for (n=first_index+1; n<=M; ++n)
    {
        const long double val = boost::math::pdf(dist, n);
        if (val > 1e-12) break;
    }

    first_index = n-1;                                                      // We start convoluting at this index, i.e., we give for granted that previous addednds where 0 as the binomial is too small (given the precision of long double)

    // Find the last value at which the binomial is larger than 1e-12 given the new value of p
    // We start from p = 1/M and n = 1 successes
    // Note that the value of n is also the last relevant index
    for (n=first_index+1; n<=M; ++n)
    {
        const long double val = boost::math::pdf(dist, n);
        binomial[n-1] = val;
        if (val < 1e-12 || n == M) break;
    }
    last_index = n;

    value = 0;
    for (n=first_index; n<=last_index; n++) { value += COP[n]*binomial[n]; }
    COP_conv[M-1] = value;

    value = 0;
    for (n=first_index; n<=last_index; n++) { value += ROP[n]*binomial[n]; }
    ROP_conv[M-1] = value;


// PRINT OUTPUT

    for(int i=0; i<M; i++) { fprintf(myOUTfile, "%.12Lf %.12Lf %.12Lf %.12Lf\n", COP[i], COP_conv[i], ROP[i], ROP_conv[i]); fflush(myOUTfile); }

    fclose(myOUTfile);


    return 0;

}
