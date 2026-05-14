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
#include "measurements.h"



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

    long int cp_mx, cp_mxy, rp_mx, rp_mxy;
    Scalars scalars;
    OrderParam CPSmax, RPSmax, Searches, TimePerBond, PivotingEvents, RigidificationEvents, OverconstrainingEvents;
    OrderParam TypeIVisited, TypeIIVisited, PivotPushes, PivotingTime, RigidificationTime, OverconstrainingTime;

    OrderParam CP_Pinf, RP_Pinf;
    LightOrderParam CP_chi, RP_chi;

    scalars.L = atoi(argv[1]);                                                  // Linear size of the lattice
    scalars.T = atoi(argv[2]);                                                  // Number of trials
    const int num = atoi(argv[3]);                                              // Run number

    basic_init(scalars, CPSmax, RPSmax, Searches, TimePerBond, PivotingEvents, RigidificationEvents, OverconstrainingEvents,
               TypeIVisited, TypeIIVisited, PivotPushes, PivotingTime, RigidificationTime, OverconstrainingTime, CP_Pinf, RP_Pinf, CP_chi, RP_chi);
    scalars.search_type = BFS;

// RELEVANT VARIABLES

    std::vector<int> RNZ(scalars.M), NZ(scalars.N);                             // Rigidity NZtree, NZ tree
    std::vector<int> bonds(scalars.M);                                          // Bonds sequence

    std::vector<std::array<int, 6>> network(scalars.N);                         // Lattice
    std::vector<std::array<int, 2>> pebble_graph(scalars.N);                    // Pebble network: pebble_graph[i] = j <=> A directed bond i->j exists; pebble_graph[i]=-1 <=> bond is inactive
    std::vector<TrialData> PerfLog(scalars.T);

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
    std::vector<int> CPdx(scalars.N);
    std::vector<int> CPdy(scalars.N);
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

    std::ostringstream LOGfname, Sfname, Searchesfname, Timefname, Eventsfname, Opsfname, EventTimesfname, PerfLogfname, CPSxyfname, RPSxyfname, WRAPfname, COPfname, ROPfname;
    LOGfname << "../res/LOG_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Execution time and generic stuff
    Sfname << "../res/Smax_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";               // Execution time and generic stuff
    Searchesfname << "../res/Searches_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    Timefname << "../res/TimePerBond_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    Eventsfname << "../res/Events_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    Opsfname << "../res/Operations_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    EventTimesfname << "../res/EventTimes_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    PerfLogfname << "../res/PERFLOG_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    CPSxyfname << "../res/CPSxy_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    RPSxyfname << "../res/RPSxy_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    WRAPfname << "../res/WRAPS_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    COPfname << "../res/COP_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";
    ROPfname << "../res/ROP_L" << scalars.L << "_T" << scalars.T << "_num" << num << ".txt";

    FILE* myLOGfile = open_file(LOGfname);
    FILE* mySfile = open_file(Sfname);
    FILE* mySearchesfile = open_file(Searchesfname);
    FILE* myTimefile = open_file(Timefname);
    FILE* myEventsfile = open_file(Eventsfname);
    FILE* myOpsfile = open_file(Opsfname);
    FILE* myEventTimesfile = open_file(EventTimesfname);
    FILE* myPerfLogfile = open_file(PerfLogfname);
    FILE* myCPSxyfile = open_file(CPSxyfname);
    FILE* myRPSxyfile = open_file(RPSxyfname);
    FILE* myWRAPfile = open_file(WRAPfname);
    FILE* myCOPfile = open_file(COPfname);
    FILE* myROPfile = open_file(ROPfname);

    fprintf(myLOGfile, "# Seed for random numers generation: %u\n", scalars.seed);
    fprintf(myLOGfile, "# L = %d, N = %d, M = %ld, T = %d\n", scalars.L, scalars.N, scalars.M, scalars.T);
    fprintf(myLOGfile, "# The execution time of each trial, expressed in seconds, follows.\n");

    fprintf(mySfile, "# Realtive size of the largest cluster for each value of m\n");
    fprintf(mySfile, "# Data are not averaged. First column: CP, second column: RP\n");
    fprintf(mySearchesfile, "# Number of pebble searches for each value of m\n");
    fprintf(mySearchesfile, "# Data are not averaged. Single column: searches\n");
    fprintf(myTimefile, "# Time needed for each bond activation\n");
    fprintf(myTimefile, "# Data are not averaged. Single column: seconds\n");
    fprintf(myEventsfile, "# Event that happened at each value of m\n");
    fprintf(myEventsfile, "# Data are not averaged. Columns: pivoting, rigidification, overconstraining\n");
    fprintf(myOpsfile, "# Number of operations for each value of m\n");
    fprintf(myOpsfile, "# Data are not averaged. Columns: type I visited, type II visited, pivot pushes\n");
    fprintf(myEventTimesfile, "# Time spent in each theorem for each value of m\n");
    fprintf(myEventTimesfile, "# Data are not averaged. Columns: pivoting, rigidification, overconstraining\n");
    fprintf(myPerfLogfile, "# Trial totals\n");
    fprintf(myPerfLogfile, "# Columns: trial time, pivoting events, rigidification events, overconstraining events, searches, type I visited, type II visited, pivot pushes, pivoting time, rigidification time, overconstraining time\n");


    fflush(myLOGfile);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// T RIGIDITY PERCOLATION TRIALS

    for(int i=0; i<scalars.T; ++i)
    {
        // To print on the screen the percentage of realizations done
        // if (i*100 % scalars.T == 0) { std::cout << "L = " << scalars.L << " num = " << num << ", i = " << i << " -> " <<  100*i/scalars.T << " %" << std::endl; }

        auto t0 = std::chrono::high_resolution_clock::now();

        init(RNZ, NZ, bonds, pebble_graph, np, P, Prc, scalars, dx, dy, CPdx, CPdy);
        single_trial(RNZ, NZ, bonds, network, pebble_graph, np, searched_bonds, marks, marks_indices, visited, visited_indices,
                     enqueued, enqueued_indices, P, Prc, ROOTS, CPSmax, RPSmax, Searches, TimePerBond, PivotingEvents, RigidificationEvents, OverconstrainingEvents,
                     TypeIVisited, TypeIIVisited, PivotPushes, PivotingTime, RigidificationTime, OverconstrainingTime, CP_Pinf, RP_Pinf, CP_chi, RP_chi, PerfLog[i],
                     scalars, dx, dy, CPdx, CPdy);

        network_reset(network, bonds, scalars);

        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<long double> dt = t1 - t0;
        PerfLog[i].trial_time = dt.count();

        find_relevant_wraps(scalars.CPwrap.X_wrap_data, scalars.CPwrap.XY_wrap_data, cp_mx, cp_mxy);
        find_relevant_wraps(scalars.RPwrap.X_wrap_data, scalars.RPwrap.XY_wrap_data, rp_mx, rp_mxy);
        fprintf(myWRAPfile, "%ld %ld %ld %ld %.4Lf %.4Lf\n", cp_mx, cp_mxy, rp_mx, rp_mxy, scalars.CP_chi_max, scalars.RP_chi_max); fflush(myWRAPfile);
        for (WrappingTriplet elem : scalars.CPwrap.XY_wrap_data) { fprintf(myCPSxyfile, "%lld\n", elem.size); fflush(myCPSxyfile); }
        for (WrappingTriplet elem : scalars.RPwrap.XY_wrap_data) { fprintf(myRPSxyfile, "%lld\n", elem.size); fflush(myRPSxyfile); }

        fprintf(myLOGfile, "%.2Lf \n", PerfLog[i].trial_time); fflush(myLOGfile);

    }
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// PRINT OUTPUT

    for(size_t i=0; i<scalars.M; i++) { fprintf(mySfile, "%.12Lf %.12LF\n", CPSmax.y[i], RPSmax.y[i]); fflush(mySfile); }
    for(size_t i=0; i<scalars.M; i++) { fprintf(mySearchesfile, "%.12Lf\n", Searches.y[i]); fflush(mySearchesfile); }
    for(size_t i=0; i<scalars.M; i++) { fprintf(myTimefile, "%.12Lf\n", TimePerBond.y[i]); fflush(myTimefile); }
    for(size_t i=0; i<scalars.M; i++) { fprintf(myEventsfile, "%.12Lf %.12Lf %.12Lf\n", PivotingEvents.y[i], RigidificationEvents.y[i], OverconstrainingEvents.y[i]); fflush(myEventsfile); }
    for(size_t i=0; i<scalars.M; i++) { fprintf(myOpsfile, "%.12Lf %.12Lf %.12Lf\n", TypeIVisited.y[i], TypeIIVisited.y[i], PivotPushes.y[i]); fflush(myOpsfile); }
    for(size_t i=0; i<scalars.M; i++) { fprintf(myEventTimesfile, "%.12Lf %.12Lf %.12Lf\n", PivotingTime.y[i], RigidificationTime.y[i], OverconstrainingTime.y[i]); fflush(myEventTimesfile); }
    for(int i=0; i<scalars.T; ++i) { fprintf(myPerfLogfile, "%.12Lf %.12Lf %.12Lf %.12Lf %.12Lf %.12Lf %.12Lf %.12Lf %.12Lf %.12Lf %.12Lf\n", PerfLog[i].trial_time, PerfLog[i].pivoting_events, PerfLog[i].rigidification_events, PerfLog[i].overconstraining_events, PerfLog[i].searches, PerfLog[i].typeI_visited, PerfLog[i].typeII_visited, PerfLog[i].pivot_pushes, PerfLog[i].pivoting_time, PerfLog[i].rigidification_time, PerfLog[i].overconstraining_time); fflush(myPerfLogfile); }
    fprintf(myCOPfile, "# %ld %ld\n", CP_Pinf.start, CP_Pinf.stop); fflush(myCOPfile);
    fprintf(myROPfile, "# %ld %ld\n", RP_Pinf.start, RP_Pinf.stop); fflush(myROPfile);
    for (size_t i=0; i<CP_Pinf.y.size(); ++i) { fprintf(myCOPfile, "%.12Lf %lld %.12Lf\n", CP_Pinf.y[i], CP_Pinf.norm[i], CP_chi.y[i]); fflush(myCOPfile); }
    for (size_t i=0; i<RP_Pinf.y.size(); ++i) { fprintf(myROPfile, "%.12Lf %lld %.12Lf\n", RP_Pinf.y[i], RP_Pinf.norm[i], RP_chi.y[i]); fflush(myROPfile); }

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    fprintf(myLOGfile, "# Max memory occupation: %ld KB\n", usage.ru_maxrss); fflush(myLOGfile);

    fclose(mySfile);
    fclose(mySearchesfile);
    fclose(myTimefile);
    fclose(myEventsfile);
    fclose(myOpsfile);
    fclose(myEventTimesfile);
    fclose(myPerfLogfile);
    fclose(myCPSxyfile);
    fclose(myRPSxyfile);
    fclose(myWRAPfile);
    fclose(myCOPfile);
    fclose(myROPfile);
    fclose(myLOGfile);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

    return 0;

}
