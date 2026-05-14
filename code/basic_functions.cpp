#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>				// Basic math functions
#include <time.h>				// Get clock time, to measure total run time
#include <iomanip>
#include <cstdlib>
#include <algorithm>

#include <bits/stdc++.h>

#include <stdlib.h>     //for using the function sleep
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <initializer_list>

#include "defs.h"
#include "basic_functions.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INITALIZATION


void basic_init(Scalars& scalars, OrderParam& CPSmax, OrderParam& RPSmax, OrderParam& Searches, OrderParam& TimePerBond,
                OrderParam& PivotingEvents, OrderParam& RigidificationEvents, OrderParam& OverconstrainingEvents,
                OrderParam& TypeIVisited, OrderParam& TypeIIVisited, OrderParam& PivotPushes,
                OrderParam& PivotingTime, OrderParam& RigidificationTime, OrderParam& OverconstrainingTime,
                OrderParam& CP_Pinf, OrderParam& RP_Pinf, LightOrderParam& CPchi, LightOrderParam& RPchi)
{

    scalars.N = scalars.L*scalars.L;
    scalars.M = 3*scalars.N;                                                    // Triangular lattice


    std::random_device rd;
    scalars.seed = rd();                                                        // TO GIVE RANDOM SEED
    // scalars.seed =1650740998;                                                // TO GIVE MANUAL SEED
    scalars.gen = std::mt19937(scalars.seed);

    CPSmax.y.resize(scalars.M, 0);
    RPSmax.y.resize(scalars.M, 0);
    Searches.y.resize(scalars.M, 0);
    TimePerBond.y.resize(scalars.M, 0);
    PivotingEvents.y.resize(scalars.M, 0);
    RigidificationEvents.y.resize(scalars.M, 0);
    OverconstrainingEvents.y.resize(scalars.M, 0);
    TypeIVisited.y.resize(scalars.M, 0);
    TypeIIVisited.y.resize(scalars.M, 0);
    PivotPushes.y.resize(scalars.M, 0);
    PivotingTime.y.resize(scalars.M, 0);
    RigidificationTime.y.resize(scalars.M, 0);
    OverconstrainingTime.y.resize(scalars.M, 0);

    CP_Pinf.start = (int)(0.28*scalars.M);
    CP_Pinf.stop = (int)(0.48*scalars.M);
    CPchi.start = CP_Pinf.start;
    CPchi.stop = CP_Pinf.stop;
    CP_Pinf.y.resize(CP_Pinf.stop + 1 - CP_Pinf.start, 0);
    CP_Pinf.norm.resize(CP_Pinf.stop + 1 - CP_Pinf.start, 0);
    CPchi.y.resize(CPchi.stop + 1 - CPchi.start, 0);

    RP_Pinf.start = (int)(0.6*scalars.M);
    RP_Pinf.stop = (int)(0.8*scalars.M);
    RPchi.start = RP_Pinf.start;
    RPchi.stop = RP_Pinf.stop;
    RP_Pinf.y.resize(RP_Pinf.stop + 1 - RP_Pinf.start, 0);
    RP_Pinf.norm.resize(RP_Pinf.stop + 1 - RP_Pinf.start, 0);
    RPchi.y.resize(RPchi.stop + 1 - RPchi.start, 0);

}




// Initialize neigh table with negative entries !! j=j+1
void network_init(std::vector<std::array<int, 6>>& network, const Scalars& scalars)
{
   for (int u=0; u<scalars.N; ++u)                                            // Each node:
   {
      network[u]={- (u + 2), -(u + scalars.L + 1), -(u + scalars.L), -u, -(u - scalars.L + 1), -(u-scalars.L+1 +1) };

      // top row
      if (u/scalars.L==scalars.L-1)
      {
         network[u][1]+=scalars.N;
         network[u][2]+=scalars.N;
      }
      // bottom row
      if (u/scalars.L==0)
      {
         network[u][4]+= -scalars.N;
         network[u][5]+= -scalars.N;
      }
      // left
      if (u%scalars.L==0)
      {
         network[u][2]+= -scalars.L;
         network[u][3]+= -scalars.L;
      }
      // right
      if ( (u+1)%scalars.L==0)
      {
         network[u][0]+= scalars.L;
         network[u][5]+= scalars.L;
      }

   }
}


void network_reset(std::vector<std::array<int, 6>>& network, std::vector<int>& bonds, const Scalars& scalars)
{
   int b, u, v, d;
   for(size_t i=0; i<scalars.M; ++i)
   {
      b = bonds[i];
      u = b/3;
      d = b%3;
      v = network[u][d] - 1;

      if( (v+1)>0 )
      {
         network[u][d]=-v-1;
         network[v][d+3]=-u-1;
      }
   }
}



void init(std::vector<int>& RNZ, std::vector<int>& NZ, std::vector<int>& bonds, std::vector<std::array<int, 2>>& pebble_graph,
          std::vector<int>& np, std::vector<int>& P, std::vector<std::set<int>>& Prc, Scalars& scalars,
          std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy, std::vector<int>& CPdx, std::vector<int>& CPdy)
{

    // Prepare datastructures for next run

    std::fill(NZ.begin(), NZ.end(), -1);                                        // Each node is its own CP root
    std::fill(RNZ.begin(), RNZ.end(), -1);                                      // Each bond is its own RP root
    std::fill(np.begin(), np.end(), 2);                                         // Each node has 2 pebbles
    std::fill(P.begin(), P.end(), 0);                                           // Each node has pivotal class 0
    std::fill(CPdx.begin(), CPdx.end(), 0);
    std::fill(CPdy.begin(), CPdy.end(), 0);

    for (size_t i=0; i<scalars.M; ++i) Prc[i] = std::set<int>();                // Rigid clusters have no pivots
    for (int i=0; i<scalars.N; ++i)
    {
        pebble_graph[i][0] = pebble_graph[i][1] = -1;                           // Pebble edges are not active
        dx[i] = std::map<int, int>();                                           // There are no displacements
        dy[i] = std::map<int, int>();
    }

    // Random order of bond activation
    std::iota(bonds.begin(), bonds.begin()+scalars.M, 0);
    std::shuffle(bonds.begin(), bonds.end(), scalars.gen);

    // Scalar variables
    scalars.n = scalars.m = 0;
    scalars.indep = scalars.red = 0;
    scalars.CPS = scalars.RPS = 1;

    scalars.x_wrap_roots.clear();
    scalars.y_wrap_roots.clear();
    scalars.xy_wrap_roots.clear();
    scalars.x_wrap = 0;
    scalars.y_wrap = 0;
    scalars.xy_wrap = 0;

    scalars.CP_chi_max = 1;
    scalars.CP_chi_num = scalars.CP_chi_den = scalars.N;
    scalars.CPwrap.X_wrap_roots.clear();
    scalars.CPwrap.Y_wrap_roots.clear();
    scalars.CPwrap.XY_wrap_roots.clear();
    scalars.CPwrap.X_wrap_data = scalars.CPwrap.Y_wrap_data = scalars.CPwrap.XY_wrap_data = std::vector<WrappingTriplet>();

    scalars.RP_chi_max = 1;
    scalars.RP_chi_num = scalars.RP_chi_den = scalars.M;
    scalars.RPwrap.X_wrap_roots.clear();
    scalars.RPwrap.Y_wrap_roots.clear();
    scalars.RPwrap.XY_wrap_roots.clear();
    scalars.RPwrap.X_wrap_data = scalars.RPwrap.Y_wrap_data = scalars.RPwrap.XY_wrap_data = std::vector<WrappingTriplet>();


}


FILE* open_file(const std::ostringstream& filenameStream) {

    std::string filename = filenameStream.str();
    FILE* file = fopen(filename.c_str(), "w");
    if (!file) {
        std::cerr << "Failed to open file " << filename << std::endl;
        exit(1);
    }
    return file;
}
