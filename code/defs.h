#ifndef DEFS_H_
#define DEFS_H_

// ALL LIBS

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <stack>
#include <cstdlib> // For RAND_MAX
#include <random>

#include <string.h>
#include <signal.h>
#include <assert.h>
#include <complex.h>

#include <ctype.h>

#define  DFS 1
#define  BFS 2

#define RIGID 1
#define FLOPPY 0
#define UNMARKED -1

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

struct WrappingTriplet {
    size_t root;
    size_t bond_num;
    unsigned long long size;
};


struct WrapState
{
    std::set<size_t> X_wrap_roots, Y_wrap_roots, XY_wrap_roots;
    std::vector<WrappingTriplet> X_wrap_data, Y_wrap_data, XY_wrap_data;
};


struct Scalars
{

	int n, m;                                                                     // Number of active nodes (n) and bonds (m)
	int indep, red;                                                               // Number of independent / redundant bonds

	int L, N, T;                                                                  // Lattice size, Number of nodes,, Number of trials
	size_t M;                                                                     // Number of Bonds

	long int CPS, RPS;                                                         // Size of the larger cluster, size of the largest rigid cluster

	int search_type;

	unsigned int seed;
	std::mt19937 gen;

	unsigned long long CP_chi_num, CP_chi_den;
	long double CP_chi_max;
	struct WrapState CPwrap;

	unsigned long long RP_chi_num, RP_chi_den;
	long double RP_chi_max;
	struct WrapState RPwrap;

	int x_wrap, y_wrap, xy_wrap;
  std::set<int> x_wrap_roots, y_wrap_roots, xy_wrap_roots;
  std::set<int> x_wrap_bond, y_wrap_bond, xy_wrap_bond;

};


struct Mixed_FIFO_LIFO {

    std::vector<int> data;
    int head = 0; 													// top of the container
    int tail = 0;                           // bottom of the container
		size_t size;

		Mixed_FIFO_LIFO(unsigned long int size) : data(size), size(size) {}

		// Insert at the bottom and increment it (both FIFO and LIFO)
		inline void push(int value) { data[tail++] = value; }


		// Remove the top element (FIFO only)
    inline void pop() { ++head; }

		// Remove the bottom element (LIFO only)
		inline void pop_back() { --tail; }


		// Access the top element (FIFO onty)
		inline int front() const { return data[head]; }

		// Access the bottom element (LIFO only)
    inline int top() const { return data[tail - 1]; }


		// Reset the container for reuse
    inline void clear() { head = tail = 0; }

		// Check if the container is empty
    inline bool empty() const { return head == tail; }

		// Expand the underlying vector
		inline void resize(unsigned long int new_size) {
			// std::cout << "I am resizing to new_size = " << new_size << std::endl;
			unsigned long int new_size_real = (800000000 > new_size) ? 800000000 : new_size;
			data.resize(new_size_real); size = new_size_real;
		}

};



struct OrderParam
{
    size_t start = 0, stop = 0;
    std::vector<long double> y;
    std::vector<unsigned long long> norm;
};


struct LightOrderParam
{
    size_t start = 0, stop = 0;
    std::vector<long double> y;
};


struct TrialData
{
    long double trial_time = 0.0;

    long double pivoting_events = 0.0;
    long double rigidification_events = 0.0;
    long double overconstraining_events = 0.0;

    long double searches = 0.0;
    long double typeI_visited = 0.0;
    long double typeII_visited = 0.0;
    long double pivot_pushes = 0.0;

    long double pivoting_time = 0.0;
    long double rigidification_time = 0.0;
    long double overconstraining_time = 0.0;
};


#endif /* DEFS_H_ */
