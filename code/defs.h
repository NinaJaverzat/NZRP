#ifndef DEFS_H_
#define DEFS_H_

// ALL LIBS

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
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



struct OrderParam { std::vector<long double> y; };


#endif /* DEFS_H_ */
