#ifndef BASIC_H_
#define BASIC_H_

#include "defs.h"

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

// INITALIZATION

void basic_init(Scalars& scalars, OrderParam& CPSmax, OrderParam& RPSmax);

void network_init(std::vector<std::array<int, 6>>& network, const Scalars& scalars);
void network_reset(std::vector<std::array<int, 6>>& network, std::vector<int>& bonds, const Scalars& scalars);


void init(std::vector<int>& RNZ, std::vector<int>& NZ, std::vector<int>& bonds, std::vector<std::array<int, 2>>& pebble_graph,
          std::vector<int>& np, std::vector<int>& P, std::vector<std::set<int>>& Prc, Scalars& scalars,
          std::vector<std::map<int, int>>& dx, std::vector<std::map<int, int>>& dy);

FILE* open_file(const std::ostringstream& filenameStream);


#endif /* BASIC_H_ */
