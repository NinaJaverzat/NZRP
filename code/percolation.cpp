#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>				    // Basic math functions
#include <chrono>           // Better clock time

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
#include "percolation.h"
#include "pebble_game.h"

static Mixed_FIFO_LIFO ToDo(1000000);
int stack[1000];


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PHASE DIAGRAM OF A SINGLE TRIAL: from empty to filled lattice

void single_trial(std::vector<int>& RNZ, std::vector<int>& NZ, const std::vector<int>& bonds, std::vector<std::array<int, 6>>& network, std::vector<std::array<int, 2>>& pebble_graph,
                  std::vector<int>& np, std::vector<int>& searched_bonds, std::vector<int>& marks, std::vector<size_t>& marks_indices,
                  std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<char>& enqueued, std::vector<size_t>& enqueued_indices,
                  std::vector<int>& P, std::vector<std::set<int>>& Prc, std::array<int, 7>& ROOTS, OrderParam& CPSmax, OrderParam& RPSmax, Scalars& scalars,
                  std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{

    int cv;

    for(int i=scalars.M-1; i>=0; --i)
    {

        // Select bond
        const int b = bonds[i];
        const int u = b/3;
        const int d = b - u*3;
        const int v = -network[u][d] - 1;

        // Roots of connectivity clusters of u and v
        const int cu = find_root_NZ(NZ, u);
        if (scalars.CPS < scalars.N) { cv = find_root_NZ(NZ, v); }
        else { cv = cu; }

        const int su = -NZ[cu];
        const int sv = -NZ[cv];

        if (cu != cv)
        {
            pivoting(RNZ, pebble_graph, network, visited, visited_indices, np, searched_bonds, P, Prc, scalars, u, v, su, sv, b, dx, dy);
            connectivity_percolation(NZ, cu, cv, scalars);
            ++scalars.indep;
        } else
        {
            const int rc_uv = same_rigid_cluster(RNZ, network, ROOTS, u, v, dx, dy);
            while (ROOTS[0]) { ROOTS[ROOTS[0]--] = -1; }

            if(rc_uv == -1)                                                     // A rigid cluster that contains both u and v does not exist
            {
                rigidification(RNZ, pebble_graph, network, searched_bonds, np, marks, marks_indices, visited, visited_indices, enqueued, enqueued_indices, P, Prc, scalars, u, v, b, dx, dy);
                ++scalars.indep;
            } else {                                                            // It exists
                overconstraining(RNZ, u, v, b, rc_uv, scalars, dx, dy);
                ++scalars.red;
            }
        }

        // Activate bond; update m
        ++scalars.m;
        network[u][d] = -network[u][d];
        network[v][d+3] = -network[v][d+3];

        update_OP(CPSmax, (long double)scalars.CPS, (long double)scalars.N, scalars.m-1);
        update_OP(RPSmax, (long double)scalars.RPS, (long double)scalars.M, scalars.m-1);

    }

    clean_visited(visited, visited_indices);                                    // unnecessary?

    if (scalars.indep - (2*scalars.N - 3) != 0) { std::cout<<"ERROR!\n"; exit(1); }
}



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// Functions for NEWMANN-ZIFF algorithm

// This function finds the root of a NZ-like tree and implements path compression
int find_root_NZ(std::vector<int>& NZ, const int node)
{

  int r, s;
  int sp=0;

  r = node;
  while ( NZ[r] >= 0 )
  {
      stack[sp++] = r;
      r = NZ[r];
  }

  while (sp)
  {
      --sp;
      s = stack[sp];
      NZ[s] = r;
  }
  return r;
}



// This function finds the root of a NZ-like tree, implements path compression
// and it further implements the update of the RP displacements
int find_root(std::vector<int>& RNZ, const int bond, std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy,
              const std::vector<std::array<int, 6>>& network, const int& v_new)
{

    int sp=0;
    int r = bond;

    while( RNZ[r] >= 0 )
    {
        stack[sp++] = r;
        r = RNZ[r];
    }

    // Compress and update the displacements along the tree
    while (sp)
    {
        --sp;
        int s=stack[sp];
        int p = RNZ[s];

       // End nodes of s -> we update their displacements wrt root r
       int u_s = s/3; int v_s = network[u_s][s%3]-1;
       if (v_s < 0) v_s = v_new;

       int u_p = p/3; // int v_p = network[u_p][p%3]-1;

       dx[u_s][r] = dx[u_s][p] + dx[u_p][r];
       dy[u_s][r] = dy[u_s][p] + dy[u_p][r];

       dx[v_s][r] = dx[v_s][p] + dx[u_p][r];
       dy[v_s][r] = dy[v_s][p] + dy[u_p][r];

       RNZ[stack[sp]] = r;

    }

    return r;

}


// Check if the rigid cluster identified by root wraps
void check_wrapping(int dx1, int dx2, int dy1, int dy2, Scalars& scalars, int root)
{

    int delta_x = abs(dx1 - dx2);
    int delta_y = abs(dy1 - dy2);

    if( delta_x>1 && delta_y<=1) // wraps horizontally
    {
        scalars.x_wrap = 1; // At least one RC percolates in x
        if (scalars.y_wrap_roots.find(root) != scalars.y_wrap_roots.end()) {scalars.xy_wrap = 1; scalars.xy_wrap_roots.insert(root); }
        scalars.x_wrap_roots.insert(root);
    }
    if( delta_x<=1 && delta_y>1) // wraps vertically
    {
        scalars.y_wrap = 1;
        if (scalars.x_wrap_roots.find(root) != scalars.x_wrap_roots.end()) {scalars.xy_wrap = 1; scalars.xy_wrap_roots.insert(root); }

        scalars.y_wrap_roots.insert(root);
    }
    if( delta_x>1 && delta_y>1)  // crossing wrap
    {
        scalars.xy_wrap=1;
        scalars.xy_wrap_roots.insert(root);
    }
}


// Merge connectivity clusters
void connectivity_percolation(std::vector<int>& NZ, const int cu, const int cv, Scalars& scalars)
{
    int large, small;

    const int su = -NZ[cu];
    const int sv = -NZ[cv];

    if (su > sv) { large =  cu; small = cv; }
    else { large = cv; small = cu; }

    NZ[large] += NZ[small];                                                     // root of "small" cluster is now "large" cluster
    NZ[small] = large;                                                          // root of "small" cluster is now "large" cluster

    if (-NZ[large] > scalars.CPS) scalars.CPS = -NZ[large];                 // If needed, update the size of the largest cluster
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// FUNCTIONS FOR RIGIDITY PERCOLATION

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// GENERIC RP FUNCTIONS

// Inreasing the pivotal class of a node might require the update of the Pivot Network
void increase_pivotal_class(std::vector<int>& RNZ, std::vector<int>& P, std::vector<std::set<int>>& Prc, const std::vector<std::array<int, 6>>& network,
                            const int node, const int b, std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    ++P[node];                                                                  // Pivotal degree grows by one

    if (P[node] > 1)
    {
        update_pivot_map(RNZ, Prc, network, node, dx, dy);        // Update pivots of rigid clusters to which u belongs, except the new one
        Prc[b].insert(node);                                             // The new one is surely pivoted by u
    }
}

// Find root of each active bond and add the node to the set of pivots
void update_pivot_map(std::vector<int>& RNZ, std::vector<std::set<int>>& Prc, const std::vector<std::array<int, 6>>& network, const int node,
                      std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    for (int d=0; d<6; ++d)
    {
        const int neigh = network[node][d] - 1;
        if (neigh+1>0)
        {
            const int bond = (d<3)*(3*node + d) + (d>=3)*(3*neigh+d-3);           // bond between node and neigh
            const int root = find_root(RNZ, bond, dx, dy, network, -1);           // We store the root of the bond that connects node to the neigh
            Prc[root].insert(node);                                               // node is a pivot of the rigid cluster identified by root
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////

// Given two input nodes, if a rigid cluster that constains both exists, the root of it is returned,
// otherwise -1 is returned
int same_rigid_cluster(std::vector<int>& RNZ, const std::vector<std::array<int, 6>>& network, std::array<int, 7>& ROOTS, const int u, const int v,
                       std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    find_all_roots(RNZ, network, ROOTS, u, dx, dy);
    return find_shared_root(RNZ, network, ROOTS, v, dx, dy);
}

// Find the root bonds of rigid clusters to which the first node (u) belongs and stores them in ROOT
void find_all_roots(std::vector<int>& RNZ, const std::vector<std::array<int, 6>>& network, std::array<int, 7>& ROOTS, const int u,
                    std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    for (int d=0; d<6; ++d)                                                     // For each neighbor of u
    {
        const int neigh = network[u][d] - 1;
        if (neigh+1>0)
        {
            const int bond = (d<3)*(3*u + d) + (d>=3)*(3*neigh+d-3);
            const int root = find_root(RNZ, bond, dx, dy, network, -1);                              // We store the root of the bond that connects u to the neighbor
            ROOTS[++ROOTS[0]] = root;
        }
    }
}

// Find the root bonds of rigid clusters to which the second node (v) belongs
// If one of them is in root, it is returned, otherwise -1 is returned
int find_shared_root(std::vector<int>& RNZ, const std::vector<std::array<int, 6>>& network, const std::array<int, 7>& ROOTS, const int v,
                     std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    for (int d=0; d<6; ++d)                                                     // For each neighbor of v
    {
        const int neigh = network[v][d] - 1;
        if (neigh+1>0)
        {
            const int bond = (d<3)*(3*v + d) + (d>=3)*(3*neigh+d-3);
            const int root = find_root(RNZ, bond, dx, dy, network, -1);                              // We compute the root of the bond between v and the neighbor

            for (int i = 1; i<=ROOTS[0]; ++i) {if (root == ROOTS[i]) { return root; }}   // If the root is also a root of u, it's returned
        }
    }
    return -1;
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PIVOTING THEOREM

// Pivoting theorem: the new bond is independent, creates a rigid cluster made of itself only
void pivoting(std::vector<int>& RNZ, std::vector<std::array<int, 2>>& pebble_graph, const std::vector<std::array<int, 6>>& network,
              std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<int>& np,
              std::vector<int>& searched_bonds, std::vector<int>& P, std::vector<std::set<int>>& Prc,
              Scalars& scalars, const int u, const int v, const int su, const int sv, const int b,
              std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    if (np[u] != 0) { cover_edge(pebble_graph, np, u, v); }                  // u has a pebble -> We use it to cover the edge
    else if (np[v] != 0) { cover_edge(pebble_graph, np, v, u); }             // u does not have pebbles, but v does -> We use it to cover the edge
    else                                                                        // A pebble must be gathered. We search for it in the smallest cluster
    {
        int small = u;
        int large = v;
        if (su > sv) { small = v; large = u; }

        pebble_game_typeI(pebble_graph, searched_bonds, np, small, visited, visited_indices, scalars);
        cover_edge(pebble_graph, np, small, large);
    }

    increase_pivotal_class(RNZ, P, Prc, network, u, b, dx, dy);
    increase_pivotal_class(RNZ, P, Prc, network, v, b, dx, dy);

    // Update displacements
    dx[u].insert({b, 0});
    dy[u].insert({b, 0});
    dx[v].insert({b, x_dist(u,v,scalars.L)});
    dy[v].insert({b, y_dist(u,v,scalars.L)});

}




/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// OVERCONSTRAINING THEOREM

// Update RNZ tree
// Checks if the activation of uv led to wrapping
void overconstraining(std::vector<int>& RNZ, const int u, const int v, const int b, const int root, Scalars& scalars,
                      std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{
    RNZ[root] += -1;
    RNZ[b] = root;

    check_wrapping(dx[u][root], dx[v][root], dy[u][root], dy[v][root], scalars, root);

    if (-RNZ[root] > scalars.RPS)  { scalars.RPS = -RNZ[root]; }
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// RIGIDIFICATION THEOREM

void rigidification(std::vector<int>& RNZ, std::vector<std::array<int, 2>>& pebble_graph, const std::vector<std::array<int, 6>>& network,
                    std::vector<int>& searched_bonds, std::vector<int>& np, std::vector<int>& marks, std::vector<size_t>& marks_indices,
                    std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<char>& enqueued, std::vector<size_t>& enqueued_indices,
                    std::vector<int>& P, std::vector<std::set<int>>& Prc, Scalars& scalars, const int u, const int v, const int new_bond,
                    std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{

    int new_root = new_bond;

    /*
    The bond is independent, hence four pebbles must be available.
    We gather them and use one of them to cover the edge.
    */
    while(np[u] < 2 && pebble_game_typeI(pebble_graph, searched_bonds, np, u, visited, visited_indices, scalars)) {;}

    visited[u] = 1; visited_indices[++visited_indices[0]] = u;
    while(np[v] < 2 && pebble_game_typeI(pebble_graph, searched_bonds, np, v, visited, visited_indices, scalars)) { visited[u] = 1; visited_indices[++visited_indices[0]] = u; }
    visited[u] = 0; visited_indices[0] = 0;

    // One of v's pebbles is used to cover the directed edge v->u
    cover_edge(pebble_graph, np, v, u);

    // At this point we are working under the ansatz that the new bond is active, is a new cluster (of current size 2) and it's activation has hence
    // increased the pivotal class of the two nodes. We will use Prc[new_bond] to BFS over the new cluster neighbors and attempt to rigidify them
    increase_pivotal_class(RNZ, P, Prc, network, u, new_bond, dx, dy);
    increase_pivotal_class(RNZ, P, Prc, network, v, new_bond, dx, dy);


    /*
    The remaining three pebbles are pinned on u and v and the mutual
    rigidity of this bond with other bonds is checked.
    The rigidi cluster is build by using a BFS over the pivot network.
    */

    // All the rigid clusters that must be coalesced are identified and coalesced
    build_new_rigid_cluster(RNZ, pebble_graph, network, searched_bonds, np, marks, marks_indices, visited, visited_indices, enqueued, enqueued_indices, P, Prc, u, v, new_root, scalars, dx, dy);

    // If the new cluster is made of uv only, the displacements are updated like in a pivoting step
    if (-RNZ[new_bond] == 1)
    {
       dx[u].insert({new_bond, 0});
       dy[u].insert({new_bond, 0});
       dx[v].insert({new_bond, x_dist(u,v,scalars.L)});
       dy[v].insert({new_bond, y_dist(u,v,scalars.L)});
    }

    if (-RNZ[new_root] > scalars.RPS) { scalars.RPS = -RNZ[new_root]; }

}



void build_new_rigid_cluster(std::vector<int>& RNZ, std::vector<std::array<int, 2>>& pebble_graph, const std::vector<std::array<int, 6>>& network,
                             std::vector<int>& searched_bonds, std::vector<int>& np, std::vector<int>& marks, std::vector<size_t>& marks_indices,
                             std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<char>& enqueued, std::vector<size_t>& enqueued_indices,
                             std::vector<int>& P, std::vector<std::set<int>>& Prc, const int u, const int v, int& new_root, Scalars& scalars,
                             std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy)
{

    int i_is_rigid, j_is_rigid;
    size_t large, small;

    int u_large, u_small, v_small;                                              // End nodes of the large and small roots, to avoid making confusion.
    int new_bond = new_root;

    // Initialization of displacement of v
    dx[v].insert({new_bond, x_dist(u,v,scalars.L)});
    dy[v].insert({new_bond, y_dist(u,v,scalars.L)});

    // Prepare the BFS over u-v rigid clusters
    // To do contains the pivots to be checked (is the stack of pivots of step 1)
    ToDo.clear();
    ToDo.push(u);
    ToDo.push(v);

    // This remains like this until the end of the BFS around u-v rigid clusters: this is how we mark u and v rigid with respect to the present bond
    marks[u] = marks[v] = RIGID;
    marks_indices[++marks_indices[0]] = u;
    marks_indices[++marks_indices[0]] = v;

    // This remains like this until the end of the BFS around u-v rigid clusters: this is how we pin the pebbles at u and v
    visited_indices[0] = 0;
    visited[u] = 1; visited_indices[++visited_indices[0]] = u;
    visited[v] = 1; visited_indices[++visited_indices[0]] = v;

    // This remains like this until the end of the BFS around u-v rigid clusters: this is how we avoid enqueueing again u and v
    enqueued_indices[0] = 0;
    enqueued[u] = 1; enqueued_indices[++enqueued_indices[0]] = u;
    enqueued[v] = 1; enqueued_indices[++enqueued_indices[0]] = v;

    // BFS
    while(!ToDo.empty())
    {
        const size_t pivot = ToDo.front(); ToDo.pop();                          // Step 1: pivot is p

        // An enqueued pivot belongs to the forming rigid cluster for sure
        if (marks[pivot] == UNMARKED)
        {
            marks[pivot] = RIGID;
            marks_indices[++marks_indices[0]] = pivot;
            pebble_graph[pivot][0] = u; pebble_graph[pivot][1] = v;             // Triangulation
        }

        for(size_t dd=0; dd<6; ++dd)                                            // Step 2 loop
        {

            const int neigh = network[pivot][dd] - 1;
            if(neigh + 1 > 0)
            {
                const int bond = (dd<3)*(3*pivot + dd) + (dd>=3)*(3*neigh+dd-3);             // This is a bond between the pivot p and one of its active neighbors
                const int root = find_root(RNZ, bond, dx, dy, network, v);                   // This is the root of such bond

                // Identify the nodes that make up the root bond
                const int i = root/3;
                int j = network[i][root - i*3]-1;
                if (j<0) { if (i==u){j=v;} else if(i==v){j=u;} }                             // This should actually not happen


                /*
                Step 3-4: CHECK IF THE ROOT OF THE CLUSTER UNDER ANALYSIS IS RIGID
                */

                i_is_rigid = j_is_rigid = FLOPPY;

                // Mark i as RIGID or FLOPPY
                if(marks[i] != UNMARKED) { i_is_rigid = marks[i]; }                     // 0 if floppy, 1 if rigid -> Only set if has been checked already
                else                                                                    // Not yet checked
                {
                    if ( np[i] != 0 )
                    {
                        marks[i] = i_is_rigid; marks_indices[++marks_indices[0]] = i;
                    } else {
                        i_is_rigid = !pebble_game_typeII (pebble_graph, searched_bonds, np, i, visited, visited_indices, marks, marks_indices, scalars);  // found = 1 -> FLOPPY and viceversa
                    }
                }

                // If i is rigid and j is UNMARKED, mark j
                if(marks[j] != UNMARKED) { j_is_rigid = marks[j]; }                     // 0 if floppy, 1 if rigid -> Only set if has been checked already
                else
                {
                    if (np[j] > 0) { marks[j] = j_is_rigid; marks_indices[++marks_indices[0]] = j;}
                    else if (i_is_rigid) {
                        j_is_rigid = !pebble_game_typeII (pebble_graph, searched_bonds, np, j, visited, visited_indices, marks, marks_indices, scalars);  // found = 1 -> FLOPPY and viceversa
                    }
                }

                /*
                  IF i is rigid and j is rigid -> the bond "root" is rigid -> all the bonds in its clusters are rigid, including pivot-neigh -> pivot and neigh are also rigid
                  NOTE: we must be careful and not update the RNZ tree if it has been updated already, i.e., if we have already redirected root (which might happen)
                */
                if (i_is_rigid && j_is_rigid && RNZ[root] != new_root && root != new_root)
                {

                    /*
                    STEP 5-8: COALESCENCE OF THE TWO RIGID CLUSTERS
                    */

                    // Preliminary: identify roots of largest and smallest rigid clusters, identify the nodes they are made of
                    u_large = i;

                    u_small = new_root/3;
                    v_small = network[u_small][new_root - u_small*3]-1;
                    if (v_small<0) { if (u_small==u){v_small=v;} else if(v_small==v){u_small=u;} }                // MAKE CHECKS. Point is: happens iff small is uv

                    if (-RNZ[root] > -RNZ[new_root]) { large =  root; small = new_root; new_root = root; }
                    else { large = new_root; small = root; u_large = u_small; u_small = i; v_small = j; }

                    // Step 5: coalescence
                    RNZ[large] += RNZ[small];
                    RNZ[small] = large;

                    // Small ceases to be a root, so we might need to remove it from the wrapping roots and replace it with large
                    if (scalars.x_wrap_roots.find(small) != scalars.x_wrap_roots.end()) {scalars.x_wrap_roots.erase(small); scalars.x_wrap_roots.insert(large);}
                    if (scalars.y_wrap_roots.find(small) != scalars.y_wrap_roots.end()) {scalars.y_wrap_roots.erase(small); scalars.y_wrap_roots.insert(large);}
                    if (scalars.xy_wrap_roots.find(small) != scalars.xy_wrap_roots.end()) {scalars.xy_wrap_roots.erase(small); scalars.xy_wrap_roots.insert(large);}

                    // Neigh belongs to the rigid cluster as one of its bonds belongs
                    if (marks[neigh] == UNMARKED)
                    {
                        marks[neigh] = RIGID;
                        marks_indices[++marks_indices[0]] = neigh;
                        pebble_graph[neigh][0] = u; pebble_graph[neigh][1] = v; // Triangulation
                    }

                    // Displacement small = (u_small, v_small) -> u_large after coalescence: Eq. (1)
                    int new_dx_us = -dx[pivot][small] + dx[pivot][large];
                    int new_dy_us = -dy[pivot][small] + dy[pivot][large];
                    int new_dx_vs = dx[v_small][small] + new_dx_us;
                    int new_dy_vs = dy[v_small][small] + new_dy_us;

                    /*
                    Check wrapping if the nodes were already pointing to the large root (this should happen only if the nodes of small are also pivots)
                    If u_small or v_small are pivots between small and large before the coalescence, we check if the coalescence of
                    large and small created a non-contractible loop passing through them.
                    */
                    if(u_small == u_large || ( dx[u_small].count(large) && dx[u_small][large] != 0 ) || ( dy[u_small].count(large) && dy[u_small][large]!=0 ))
                    {
                        check_wrapping(new_dx_us, dx[u_small][large], new_dy_us, dy[u_small][large], scalars, large);
                    }
                    if(v_small == u_large || ( dx[v_small].count(large) && dx[v_small][large] != 0 ) || ( dy[v_small].count(large) && dy[v_small][large]!=0 ))
                    {
                        check_wrapping(new_dx_vs, dx[v_small][large], new_dy_vs, dy[v_small][large], scalars, large);
                    }

                    // Store the new displacements
                    dx[u_small][large] = new_dx_us;
                    dy[u_small][large] = new_dy_us;

                    dx[v_small][large] = new_dx_vs;
                    dy[v_small][large] = new_dy_vs;


                    // Check wrapping from pivot
                    check_wrapping(dx[pivot][large], dx[pivot][small]+dx[u_small][large], dy[pivot][large], dy[pivot][small]+dy[u_small][large], scalars, large);

                    /*
                    Step 6-8: UPDATE OF THE PIVOT MAP
                    */

                    // Step 6
                    --P[pivot];
                    Prc[small].erase(pivot);
                    if (P[pivot] == 1) Prc[large].erase(pivot);

                    // Step 7: piv is p'
                    for (int piv : Prc[small])
                    {

                        // Displacement piv -> u_large after coalescence: Eq. (2)
                        int new_dx = dx[piv][small] + dx[u_small][large];
                        int new_dy = dy[piv][small] + dy[u_small][large];

                        if (Prc[large].find(piv) != Prc[large].end())           // If condition of step 7
                        {
                            --P[piv];
                            if(P[piv]==1) Prc[large].erase(piv);
                            check_wrapping(new_dx, dx[piv][large], new_dy, dy[piv][large], scalars, large);        // Check
                        } else {                                                 // Still step 7
                            Prc[large].insert(piv);
                            if (enqueued[piv] == 0)                              // Step 8
                            {
                                enqueued[piv] = 1;
                                enqueued_indices[++enqueued_indices[0]] = piv;
                                ToDo.push(piv);
                            }
                        }

                        // Store the new displacements
                        dx[piv][large] = new_dx;
                        dy[piv][large] = new_dy;
                    }
                    Prc[small] = std::set<int>();                               // small is not a root anymore
                }
            }
        }
    }

    // Reinitialization of marks
    clean_marks(marks, marks_indices);
    clean_visited(visited, visited_indices);                                    // Visited is actually supposed to be clean when we enter here
    clean_visited(enqueued, enqueued_indices);

}















//
