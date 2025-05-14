#ifndef PERC_H_
#define PERC_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PERCOLATION FUNCTIONS

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


// PHASE DIAGRAM OF A SINGLE TRIAL

void single_trial(std::vector<int>& RNZ, std::vector<int>& NZ, const std::vector<int>& bonds, std::vector<std::array<int, 6>>& network, std::vector<std::array<int, 2>>& pebble_graph,
                  std::vector<int>& np, std::vector<int>& searched_bonds, std::vector<int>& marks, std::vector<size_t>& marks_indices,
                  std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<char>& enqueued, std::vector<size_t>& enqueued_indices,
                  std::vector<int>& P, std::vector<std::set<int>>& Prc, std::array<int, 7>& ROOTS, OrderParam& CPSmax, OrderParam& RPSmax, Scalars& scalars,
                  std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// Functions for NEWMANN-ZIFF algorithm


int find_root_NZ(std::vector<int>& NZ, const int node);

int find_root(std::vector<int>& RNZ, const int bond, std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy,
              const std::vector<std::array<int, 6>>& network, const int& v_new);

void check_wrapping(int dx1, int dx2, int dy1, int dy2, Scalars& scalars, int root);

void connectivity_percolation(std::vector<int>& NZ, const int cu, const int cv, Scalars& scalars);

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// GENERIC RP FUNCTIONS

void increase_pivotal_class(std::vector<int>& RNZ, std::vector<int>& P, std::vector<std::set<int>>& Prc, const std::vector<std::array<int, 6>>& network,
                            const int node, const int b, std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

void update_pivot_map(std::vector<int>& RNZ, std::vector<std::set<int>>& Prc, const std::vector<std::array<int, 6>>& network, const int node,
                      std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

int same_rigid_cluster(std::vector<int>& RNZ, const  std::vector<std::array<int, 6>>& network, std::array<int, 7>& ROOTS, const int u, const int v,
                       std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

void find_all_roots(std::vector<int>& RNZ, const std::vector<std::array<int, 6>>& network, std::array<int, 7>& ROOTS, const int u,
                    std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

int find_shared_root(std::vector<int>& RNZ, const std::vector<std::array<int, 6>>& network, const std::array<int, 7>& ROOTS, const int v,
                     std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// PIVOTING THEOREM

void pivoting(std::vector<int>& RNZ, std::vector<std::array<int, 2>>& pebble_graph, const std::vector<std::array<int, 6>>& network,
              std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<int>& np,
              std::vector<int>& searched_bonds, std::vector<int>& P, std::vector<std::set<int>>& Prc,
              Scalars& scalars, const int u, const int v, const int su, const int sv, const int b,
              std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// OVERCONSTRAINING THEOREM

void overconstraining(std::vector<int>& RNZ, const int u, const int v, const int b, const int root, Scalars& scalars,
                      std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// RIGIDIFICATION THEOREM

void rigidification(std::vector<int>& RNZ, std::vector<std::array<int, 2>>& pebble_graph, const std::vector<std::array<int, 6>>& network,
                    std::vector<int>& searched_bonds, std::vector<int>& np, std::vector<int>& marks, std::vector<size_t>& marks_indices,
                    std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<char>& enqueued, std::vector<size_t>& enqueued_indices,
                    std::vector<int>& P, std::vector<std::set<int>>& Prc, Scalars& scalars, const int u, const int v, const int new_bond,
                    std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);


void build_new_rigid_cluster(std::vector<int>& RNZ, std::vector<std::array<int, 2>>& pebble_graph, const std::vector<std::array<int, 6>>& network,
                             std::vector<int>& searched_bonds, std::vector<int>& np, std::vector<int>& marks, std::vector<size_t>& marks_indices,
                             std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<char>& enqueued, std::vector<size_t>& enqueued_indices,
                             std::vector<int>& P, std::vector<std::set<int>>& Prc, const int u, const int v, int& new_root, Scalars& scalars,
                             std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INLINE FUNCTIONS

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INLINE FUNCTIONS

// HANDLING OF THE VISITED AND MARKS ARRAYS

inline void clean_marks(std::vector<int>& marks, std::vector<size_t>& marks_indices)
{
    while(marks_indices[0] >= 4){
      marks[marks_indices[marks_indices[0]--]] = UNMARKED;
      marks[marks_indices[marks_indices[0]--]] = UNMARKED;
      marks[marks_indices[marks_indices[0]--]] = UNMARKED;
      marks[marks_indices[marks_indices[0]--]] = UNMARKED;
    }
    while(marks_indices[0]){ marks[marks_indices[marks_indices[0]--]] = UNMARKED; }
}


inline void clean_visited(std::vector<char>& visited, std::vector<size_t>& visited_indices)
{
    while(visited_indices[0] >= 4){
      visited[visited_indices[visited_indices[0]--]] = 0;
      visited[visited_indices[visited_indices[0]--]] = 0;
      visited[visited_indices[visited_indices[0]--]] = 0;
      visited[visited_indices[visited_indices[0]--]] = 0;
    }
    while(visited_indices[0]){visited[visited_indices[visited_indices[0]--]] = 0;}
}

inline void clean_part_of_visited(std::vector<char>& visited, std::vector<size_t>& visited_indices)
{
  while(visited_indices[0] >= 6){
    visited[visited_indices[visited_indices[0]--]] = 0;
    visited[visited_indices[visited_indices[0]--]] = 0;
    visited[visited_indices[visited_indices[0]--]] = 0;
    visited[visited_indices[visited_indices[0]--]] = 0;
  }
  while(visited_indices[0] > 2){visited[visited_indices[visited_indices[0]--]] = 0;}
}

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// Flips the sign of the displacement if the bond between u and v crosses the boundary:
// - Returns -1 if the boundary is crossed from left to right
// - Returns 1 if the boundary is crossed from right to left
inline int dist_sign(const int x){ if (abs(x)<=1) return x; else return ((x < 0) - (x > 0)); }

// Computes the x and y algebraic displacements from v to u, in order to update the Machta displacements.
// u and v are always neighbours, ie are distant by one lattice spacing so that these functions return +/- 1 or 0.
// The function dist_sign is used to account for periodic boundary crossing. 
inline int x_dist(const int u, const int v, const int& L) { return dist_sign( u%L - v%L );}
inline int y_dist(const int u, const int v, const int& L) { return dist_sign( u/L - v/L );}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// Computes the number of pebbles from the pebble graph
// inline int np(const std::array<int, 2>& pebble_bonds) {
//
//   if (pebble_bonds[1] != -1) return 0;
//   if (pebble_bonds[0] == -1) return 2;
//   return 1;
// }

// Update the pebble graph and the number of pebbles
inline void cover_edge( std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& np, const int from, const int to)
{
    pebble_graph[from][2 - np[from]] = to;
    --np[from];
}


inline void update_OP(OrderParam& OP, const long double val, const long double norm, const int ind)
{
    OP.y[ind] += val / norm;
}



#endif /* PERC_H_ */
