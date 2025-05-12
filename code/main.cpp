/*
Bond-diluted rigidity percolation with pivot-based NZ approach.

1) Activate bonds in random order;
2) Study connectivity percolation using the NZ algorithm;
3) If:
A) Different connectivity cluster -> Pivoting -> New rigid cluster made of the new bond only;
B) Same connectivity cluster + same rigid cluster -> Overconstaining -> Redundant bond;
C) Same connectivity cluster + different rigid cluster -> Rigidification (Independent bond) ->

--- Rigidification is performed by
--- iterating over pivots node of the forming rigid cluster. Initially the pivots are just the
--- two nodes connected by the new bond, which are guaranteed to belong to other rigid clusters as well.
--- Trough the pivots, the root bonds of these clusters are identified and checked. If they
--- are mutually rigid with respect to the new bond, the Newmann-Ziff tree and the list of pivots is updated
--- accordingly. The new pivots, then, will be cheked until all the rigidifying clusters are identified
--- and merged onto the new one. Crucially, when we merge two clusters. we merge the small one in the large one,
--- as it (typically) has much less pivots, which makes the merging faster.

--- See: LINK TO GIT, LINK TO PAPER

*/

#include <cstdio>
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <initializer_list>

#include "defs.h"
#include "basic_functions.h"
#include "percolation.h"



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// MAIN FUNCTION

int main(int argc, char* argv[]) {

    if (argc!=4) {
      std::cerr << "Usage: " << argv[0] << " <LATTICE LINEAR SIZE>  <NUMBER OF SAMPLES> <NUMBER OF THE RUN> " << std::endl;
      return 1;
    }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// BASIC VARIABLES

    Scalars scalars;
    OrderParam CPSmax, RPSmax;                                                  // Relative size of largest CP and RP cluster respectively

    scalars.L = atoi(argv[1]);                                                  // Linear size of the lattice
    scalars.T = atoi(argv[2]);                                                  // Number of trials
    const int num = atoi(argv[3]);                                              // Run number

    basic_init(scalars, CPSmax, RPSmax);
    scalars.search_type = BFS;

// RELEVANT VARIABLES

    std::vector<int> RNZ(scalars.M), NZ(scalars.N);                             // Rigidity NZtree, NZ tree
    std::vector<int> bonds(scalars.M);                                          // Bonds sequence

    std::vector<std::array<int, 6>> network(scalars.N);                         // Lattice
    std::vector<std::array<int, 2>> pebble_graph(scalars.N);                    // Pebble network: pebble_graph[i] = j <=> A directed bond i->j exists; pebble_graph[i]=-1 <=> bond is inactive

    // Path reconstructor for BFS and DFS
    std::vector<int> searched_bonds(scalars.N);

    // Number of pebbles per node
    std::vector<int> np(scalars.N);

    // Node state marker: RIGID, FLOPPY, UNMARKED
    std::vector<int> marks(scalars.N, UNMARKED);                                // Marks to put on nodes duing game type II
    std::vector<size_t> marks_indices(scalars.N+1);                             // Filled from the head, the first element counts where we are writing

    // Visited nodes marker: visited or not
    std::vector<char> visited(scalars.N, 0);
    std::vector<size_t> visited_indices(scalars.N+1);                           // Filled from the head, the first element counts where we are writing

    // Enqueued nodes marker: enqueued or not
    std::vector<char> enqueued(scalars.N, 0);
    std::vector<size_t> enqueued_indices(scalars.N+1);                          // Filled from the head, the first element counts where we are writing

    // Pivotal class of nodes and pivots of rigid clusters
    std::vector<int> P(scalars.N, 0);
    std::vector<std::set<int>> Prc(scalars.M);                                  // This is the Pivot Network

    // Wrapping
    std::vector<std::map<int, int>> dx(scalars.N);			                      	// dx[r1][r2] = displacement from bond r1 to the first node of bond r2
    std::vector<std::map<int, int>> dy(scalars.N);

    // Array containing root bonds of a node
    std::array<int, 7> ROOTS;

    ROOTS[0] = 0;
    marks_indices[0] = 0;
    network_init(network, scalars);



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// OUTPUT FILES

    std::ostringstream LOGfname, Sfname;
    LOGfname << "../res/LOG_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Execution time and generic stuff
    Sfname << "../res/Smax_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Execution time and generic stuff

    FILE* myLOGfile = open_file(LOGfname);
    FILE* mySfile = open_file(Sfname);

    fprintf(myLOGfile, "# Seed for random numers generation: %u\n", scalars.seed);
    fprintf(myLOGfile, "# L = %d, N = %d, M = %ld, T = %d\n", scalars.L, scalars.N, scalars.M, scalars.T);
    fprintf(myLOGfile, "# The execution time of each trial, expressed in seconds, follows.\n");

    fprintf(mySfile, "# Realtive size of the largest cluster for each value of m\n");
    fprintf(mySfile, "# Data are not averaged. First column: CP, second column: RP\n");


    fflush(myLOGfile);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// T RIGIDITY PERCOLATION TRIALS

    for(int i=0; i<scalars.T; ++i)
    {
        // To print on the screen the percentage of realizations done
        // if (i*100 % scalars.T == 0) { std::cout << "L = " << scalars.L << " num = " << num << ", i = " << i << " -> " <<  100*i/scalars.T << " %" << std::endl; }

        auto t0 = std::chrono::high_resolution_clock::now();

        init(RNZ, NZ, bonds, pebble_graph, np, P, Prc, scalars, dx, dy);
        single_trial(RNZ, NZ, bonds, network, pebble_graph, np, searched_bonds, marks, marks_indices, visited, visited_indices,
                     enqueued, enqueued_indices, P, Prc, ROOTS, CPSmax, RPSmax, scalars, dx, dy);

        network_reset(network, bonds, scalars);

        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = t1 - t0;

        fprintf(myLOGfile, "%.2f \n", dt.count()); fflush(myLOGfile);

    }
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// PRINT OUTPUT

    for(size_t i=0; i<scalars.M; i++) { fprintf(mySfile, "%.12Lf %.12LF\n", CPSmax.y[i], RPSmax.y[i]); fflush(mySfile); }

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    fprintf(myLOGfile, "# Max memory occupation: %ld KB\n", usage.ru_maxrss); fflush(myLOGfile);

    fclose(mySfile);
    fclose(myLOGfile);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    return 0;

}
