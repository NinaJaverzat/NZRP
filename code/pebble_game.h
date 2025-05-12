#ifndef PG_H_
#define PG_H_

#include "defs.h"


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// PEBBLE GAME FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// TYPE 1 GAME

bool pebble_game_typeI ( std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np, const int start,
                          std::vector<char>& visited, std::vector<size_t>& visited_indices, Scalars& scalars);


bool find_path_typeI_DFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                         std::vector<char>& visited, std::vector<size_t>& visited_indices, const int start, int *endpoint);

bool find_path_typeI_BFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                        std::vector<char>& visited, std::vector<size_t>& visited_indices, const int start, int* endpoint);

void reverse_path(std::vector<int>& searched_bonds, std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& np, const int start, const int endpoint);

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// TYPE 2 GAME

bool pebble_game_typeII (std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds,
                         std::vector<int>& np, const int start, std::vector<char>& visited, std::vector<size_t>& visited_indices,
                         std::vector<int>& marks, std::vector<size_t>& marks_indices, Scalars& scalars);

bool find_path_typeII_DFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                         std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<int>& marks, const int start, int* endpoint);


bool find_path_typeII_BFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                         std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<int>& marks, const int start, int* endpoint);


void triangluate_and_mark_rigid(std::vector<std::array<int, 2>>& pebble_graph, std::vector<char>& visited, std::vector<size_t>& visited_indices,
                                std::vector<int>& marks, std::vector<size_t>& marks_indices, const int u, const int v);


inline void mark_floppy(std::vector<int>& searched_bonds, std::vector<int>& marks, std::vector<size_t>& marks_indices, const int start, const int endpoint)
{
    int node = endpoint;
    marks[node] = FLOPPY; marks_indices[++marks_indices[0]] = node;
    while(node != start) { node = searched_bonds[node]; marks[node] = FLOPPY; marks_indices[++marks_indices[0]] = node; }
}



#endif /* PG_H_ */
