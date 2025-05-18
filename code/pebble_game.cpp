#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>				// Basic math functions
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
#include "pebble_game.h"
#include "percolation.h"

static Mixed_FIFO_LIFO searcher(10000000);

/////////////////////////////////////////////////////////////////////////////////
// In this block of functions we play the pebble game.
/////////////////////////////////////////////////////////////////////////////////


/*
TYPE 1 pebble game:
- We dont "search" a pebble, we know we can find it, we simply go and get it.
- The path to the pebble is reversed.
- If the pebble is not found, something is wrong.
- In particular, as the pebble is always found, we never triangulate.
- Marks do not exist and indeed we do not mark nodes as rigid or floppy.
*/


bool pebble_game_typeI (std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np, const int start,
                        std::vector<char>& visited, std::vector<size_t>& visited_indices, Scalars& scalars)
{
    bool found = 0;
    int endpoint;

    if (scalars.search_type == DFS) { found = find_path_typeI_DFS(pebble_graph, searched_bonds, np, visited, visited_indices, start, &endpoint); }
    else { found = find_path_typeI_BFS(pebble_graph, searched_bonds, np, visited, visited_indices, start, &endpoint); }

    reverse_path(searched_bonds, pebble_graph, np, start, endpoint);
    clean_visited(visited, visited_indices);

    return found;
}

/////////////////////////////////////////////////////////////////////////////////
// Type 1 DFS

bool find_path_typeI_DFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                         std::vector<char>& visited, std::vector<size_t>& visited_indices, const int start, int *endpoint)
{

    searcher.clear();                                                           // Here we use the container like a stack. It starts empty
    searcher.push(start);                                                       // Insert at the end

    visited[start] = 1;                                                         // The starting node is a visited one
    visited_indices[++visited_indices[0]] = start;

    while (!searcher.empty())                                                   // Until we are forced to retreat all the way back to the first vertex...
    {
        const int node = searcher.top();                                        // Extract from the tail (like a stack)

        for (size_t i=0; i<2; ++i)                                              // Check each pebble neighbor of the current node
        {
            const int next_node = pebble_graph[node][i];
            if (next_node == -1) break;

            if (!visited[next_node])                                            // If the node has not been checked yet in this search, we check it and take note of the check
            {
                visited[next_node] = 1;
                visited_indices[++visited_indices[0]] = next_node;
                searched_bonds[next_node] = node;

                if (np[next_node] > 0) { *endpoint = next_node; return 1; }     // If the node has a pebble we leave

                searcher.push(next_node);
                if ( (size_t)searcher.tail == searcher.size ) searcher.resize(searcher.size*1.5);
                break;
            }
        }

        if (node == searcher.top()) searcher.pop_back();                        // We retreat if we have reached the deepest unvisited node in the branch (DFS)

    }
    return 0;
}


/////////////////////////////////////////////////////////////////////////////////
// Type 1 BFS

bool find_path_typeI_BFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                         std::vector<char>& visited, std::vector<size_t>& visited_indices, const int start, int* endpoint)
{

    searcher.clear();                                                           // Here we use the container like a queue. It start empty
    searcher.push(start);                                                       // Insert at the end

    visited[start] = 1;                                                         // The starting node is a visited one
    visited_indices[++visited_indices[0]] = start;

    while (!searcher.empty())                                                   // Until we are forced to retreat all the way back to the first vertex...
    {
        const int node = searcher.front();                                      // Extract from the head (like a queue) and retreat (BFS)
        searcher.pop();

        for (size_t i=0; i<2; ++i)                                              // Check each pebble neighbor of the current node
        {
            const int next_node = pebble_graph[node][i];

            if (next_node == -1) break;

            if (!visited[next_node])                                            // If the node has not been checked yet in this search, we check it and take note of the check
            {
                visited[next_node] = 1;
                visited_indices[++visited_indices[0]] = next_node;
                searched_bonds[next_node] = node;

                if (np[next_node] > 0) { *endpoint = next_node; return 1; }     // If the node has a pebble we leave

                searcher.push(next_node);
                if ( (size_t)searcher.tail == searcher.size ) searcher.resize(searcher.size*1.5);
            }
        }
    }
    return 0;
}


void reverse_path(std::vector<int>& searched_bonds, std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& np, const int start, const int endpoint)
{

    int node = endpoint;
    --np[node];

    while(node != start)
    {
        const int next_node = searched_bonds[node];                             // Backward tracing

        if (pebble_graph[node][0] == -1) { pebble_graph[node][0] = next_node; } // Put in cell 0 only if it is empty
        else { pebble_graph[node][1] = next_node; }                             // Otherwise put it in cell 1

        // If the node that has to be forgotten is the first, we shift the one to be rememebered in the first position
        if (pebble_graph[next_node][0] == node) pebble_graph[next_node][0] = pebble_graph[next_node][1];

        // If the node has two directed edges -> the second of them is erased
        if (pebble_graph[next_node][1] != -1) { pebble_graph[next_node][1] = -1; }
        else { pebble_graph[next_node][0] = -1; }                              // Otherwise it has one and that is the one that gets erased

        node = next_node;
    }
    ++np[node];
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

/*
TYPE 2 pebble game:
- This type of game is only played during rigidification.
- We search a pebble but do not need to find it, we typically want to discover if it is available or not.
- The path to the pebble is NOT reversed as we do not move pebbles even if they are found.
- Marks do exists: if marks[node] == FLOPPY/RIGID then node is floppy/rigid, otherwise it is UNMARKED.
- Thanks to this, we only add to visited nodes that are UNMARKED.
- If the pebble is not found, we triangulate all of visited nodes over a base and mark them RIGID.
- The base is always made of the the two bonds that are connected by the newly activated bond that triggered the rigidification.
- If the pebble is found (or a FLOPPY node is encountered), the path to the pebble (or the the FLOPPY node) is marked as FLOPPY.
*/


bool pebble_game_typeII (std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds,
                         std::vector<int>& np, const int start, std::vector<char>& visited, std::vector<size_t>& visited_indices,
                         std::vector<int>& marks, std::vector<size_t>& marks_indices, Scalars& scalars)
{
    int u, v, endpoint;
    bool found = 0;

    u = visited_indices[1];
    v = visited_indices[2];

    if (scalars.search_type == DFS) { found = find_path_typeII_DFS(pebble_graph, searched_bonds, np, visited, visited_indices, marks, start, &endpoint); }
    else { found = find_path_typeII_BFS(pebble_graph, searched_bonds, np, visited, visited_indices, marks,  start, &endpoint); }

    if (found)
    {
        mark_floppy (searched_bonds, marks, marks_indices, start, endpoint);
        clean_part_of_visited(visited, visited_indices);

        return 1;
    }

    // Basic
    triangluate_and_mark_rigid (pebble_graph, visited, visited_indices, marks, marks_indices, u, v);

    return 0;
}


/////////////////////////////////////////////////////////////////////////////////
// Type 2 DFS


bool find_path_typeII_DFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                          std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<int>& marks, const int start, int* endpoint)
{

    searcher.clear();                                                           // Here we use the container like a stack. It start being empty
    searcher.push(start);                                                       // Insert at the end

    visited[start] = 1;
    visited_indices[++visited_indices[0]] = start;

    while (!searcher.empty())                                                   // Until we are forced to retreat all the way back to the first vertex...
    {
        const int node = searcher.top();                                        // Extract from the tail (like a stack)

        for (size_t i=0; i<2; ++i)                                              // Check each neighbor of the current node
        {
            const int next_node = pebble_graph[node][i];
            if (next_node == -1) break;

            // (i) We know already that it is floppy -> We will find a pebble for sure -> the whole path to it is floppy
            if (marks[next_node] == FLOPPY) { *endpoint = node; return 1; }

            // (ii):
            // - We do not enter if it is rigid
            // - We do not enter if it has been visited already
            // - We do not insert in visited a node that is already maked rigid (and hence necessarily already triangulated)
            if(!visited[next_node] && marks[next_node] != RIGID)
            {
                visited[next_node] = 1;
                visited_indices[++visited_indices[0]] = next_node;
                searched_bonds[next_node] = node;

                if (np[next_node] > 0) { *endpoint = next_node; return 1; }

                searcher.push(next_node);
                if ( (size_t)searcher.tail == searcher.size ) searcher.resize(searcher.size*1.5);
                break;
            }
        }

        if (node == searcher.top()) searcher.pop_back();                        // We retreat if we have reached the deepest unvisited node in the branch (DFS)

    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////
// Type 2 BFS


bool find_path_typeII_BFS(const std::vector<std::array<int, 2>>& pebble_graph, std::vector<int>& searched_bonds, std::vector<int>& np,
                          std::vector<char>& visited, std::vector<size_t>& visited_indices, std::vector<int>& marks, const int start, int* endpoint)
{

    searcher.clear();                                                           // Here we use the container like a queue. It start being empty
    searcher.push(start);                                                       // Insert at the end

    visited[start] = 1;
    visited_indices[++visited_indices[0]] = start;

    while (!searcher.empty())                                                    // Until we are forced to retreat all the way back to the first vertex...
    {
        const int node = searcher.front();                                      // Extract from the head (like a queue) and retreat (BFS)
        searcher.pop();

        for (size_t i=0; i<2; ++i)                                              // Check each neighbor of the current node
        {
            const int next_node = pebble_graph[node][i];
            if (next_node == -1) break;

            // (i) We know already that it is floppy -> We will find a pebble for sure -> the whole path to it is floppy
            if (marks[next_node] == FLOPPY) { *endpoint = node; return 1; }

            // (ii):
            // - We do not enter if it is rigid
            // - We do not enter if it has been visited already
            // - We do not insert in visited a node that is already maked rigid (and hence necessarily already triangulated)
            if(!visited[next_node] && marks[next_node] != RIGID)
            {
                visited[next_node] = 1;
                visited_indices[++visited_indices[0]] = next_node;
                searched_bonds[next_node] = node;

                if (np[next_node] > 0) { *endpoint = next_node; return 1; }

                searcher.push(next_node);
                if ( (size_t)searcher.tail == searcher.size ) searcher.resize(searcher.size*1.5);
            }
        }
    }
    return 0;
}

// Triangulation: the pebble edges of the visited nodes, that are surely both active, are redirected
// Visited nodes are marked rigid
void triangluate_and_mark_rigid(std::vector<std::array<int, 2>>& pebble_graph, std::vector<char>& visited, std::vector<size_t>& visited_indices,
                                std::vector<int>& marks, std::vector<size_t>& marks_indices, const int u, const int v)
{

    for(size_t i=3; i<=visited_indices[0]; ++i)
    {
        const size_t node = visited_indices[i];

        marks[node] = RIGID;
        marks_indices[++marks_indices[0]] = node;

        pebble_graph[node][0] = u;
        pebble_graph[node][1] = v;
        visited[node] = 0;                                           // All visited nodes gets cleared inside visited, except u-v
    }
    visited_indices[0] = 2;                                          // visited_indices are reset so to keep u-v rigid and visited

}
