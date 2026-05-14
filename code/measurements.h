#ifndef MEAS_H_
#define MEAS_H_

#include "defs.h"
#include <iostream>
#include <optional>

void update_chi(const WrapState& wrap, const int& large, const int& small, const unsigned long long& large_size,
                const unsigned long long& small_size, unsigned long long& chi_num, unsigned long long& chi_den);

void wrapping(const std::vector<int>& ClusterTree, const int& dx1, const int& dx2, const int& dy1, const int& dy2, const size_t& root,
              const size_t& bond_num, WrapState& wrap, unsigned long long& chi_num, unsigned long long& chi_den, const unsigned long long tot);

void measurements(const std::vector<int>& NZ, const std::vector<int>& RNZ, OrderParam& CP_Pinf, OrderParam& RP_Pinf,
                  LightOrderParam& CP_chi, LightOrderParam& RP_chi, Scalars& scalars);

void check_master_equation(const std::vector<int>& ClusterTree, const WrapState& wrap, const unsigned long long& chi_den, const unsigned long long tot);

void find_relevant_wraps(const std::vector<WrappingTriplet>& x_wrap, const std::vector<WrappingTriplet>& xy_wrap, long int& x_first_bond, long int& xy_first_bond);

inline std::optional<size_t> find_index_of_wrapping_root(const std::vector<WrappingTriplet>& wrap_data, const size_t& root)
{
    for (size_t i=0; i<wrap_data.size(); ++i) { if (wrap_data[i].root == root) return i; }
    return std::nullopt;
}

inline void new_observed_wrap(std::set<size_t>& wrap_roots, std::vector<WrappingTriplet>& wrap_data,
                        const size_t& root, const size_t& bond_num, const unsigned long long& size)
{
    if (wrap_roots.empty() || wrap_roots.find(root) == wrap_roots.end())
    {
        wrap_roots.insert(root);
        wrap_data.emplace_back(WrappingTriplet{.root = root, .bond_num = bond_num, .size = size});
    }
}

inline void check_cross_wrap(std::set<size_t>& wrap_roots_to_check, std::set<size_t>& wrap_roots_to_udpate, std::vector<WrappingTriplet>& wrap_data,
                        const size_t& root, const size_t& bond_num, const unsigned long long& size)
{
  if (wrap_roots_to_check.find(root) != wrap_roots_to_check.end())
  {
      auto root_index = find_index_of_wrapping_root(wrap_data, root);
      if (!root_index.has_value()) { wrap_data.emplace_back(WrappingTriplet{.root = root, .bond_num = bond_num, .size = size}); wrap_roots_to_udpate.insert(root); }
  }
}

inline void update_wrap(std::set<size_t>& wrap_roots, std::vector<WrappingTriplet>& wrap_data, const size_t& large, const size_t& small)
{
    if (wrap_roots.find(small) != wrap_roots.end())
    {
        auto root_index = find_index_of_wrapping_root(wrap_data, small);
        if (root_index.has_value()) { wrap_data[root_index.value()].root = large; }
        else { std::cout << "PROBLEM! A root was in wrap_roots but it was not in wrap_data!" <<std::endl; }
        wrap_roots.erase(small); wrap_roots.insert(large);
    }
}

inline void update_OP_light(LightOrderParam& OP, const long double& val, const long double& norm, const size_t& ind)
{
    OP.y[ind] += val / norm;
}

#endif /* MEAS_H_ */
