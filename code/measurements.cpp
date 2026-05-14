#include <iostream>
#include <vector>
#include <stack>
#include <fstream>
#include <math.h>

#include <iomanip>
#include <cstdlib>
#include <algorithm>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#include <initializer_list>

#include "defs.h"
#include "percolation.h"
#include "measurements.h"
#include "basic_functions.h"

void update_chi(const WrapState& wrap, const int& large, const int& small, const unsigned long long& large_size,
                const unsigned long long& small_size, unsigned long long& chi_num, unsigned long long& chi_den)
{
    if ( (wrap.XY_wrap_roots.empty() || wrap.XY_wrap_roots.find(small) == wrap.XY_wrap_roots.end()) &&
         (wrap.XY_wrap_roots.empty() || wrap.XY_wrap_roots.find(large) == wrap.XY_wrap_roots.end()) )
    {
        chi_num += 2*large_size*small_size;
    }
    else if ( (wrap.XY_wrap_roots.empty() || wrap.XY_wrap_roots.find(small) == wrap.XY_wrap_roots.end()) &&
              (wrap.XY_wrap_roots.find(large) != wrap.XY_wrap_roots.end()) )
    {
        chi_num -= (small_size * small_size);
        chi_den -= (small_size);
    }
    else if ( (wrap.XY_wrap_roots.find(small) != wrap.XY_wrap_roots.end()) &&
              (wrap.XY_wrap_roots.empty() || wrap.XY_wrap_roots.find(large) == wrap.XY_wrap_roots.end()) )
    {
        chi_num -= large_size * large_size;
        chi_den -= large_size;
    }
}


void wrapping(const std::vector<int>& ClusterTree, const int& dx1, const int& dx2, const int& dy1, const int& dy2, const size_t& root,
              const size_t& bond_num, WrapState& wrap, unsigned long long& chi_num, unsigned long long& chi_den, const unsigned long long tot)
{
    const int delta_x = abs(dx1 - dx2);
    const int delta_y = abs(dy1 - dy2);

    bool remove_from_chi = 1;
    if (wrap.XY_wrap_roots.find(root) != wrap.XY_wrap_roots.end()) remove_from_chi = 0;

    if( delta_x>1 && delta_y<=1)
    {
        new_observed_wrap(wrap.X_wrap_roots, wrap.X_wrap_data, root, bond_num, -ClusterTree[root]);
        check_cross_wrap(wrap.Y_wrap_roots, wrap.XY_wrap_roots, wrap.XY_wrap_data, root, bond_num, -ClusterTree[root]);
    }

    if( delta_x<=1 && delta_y>1)
    {
        new_observed_wrap(wrap.Y_wrap_roots, wrap.Y_wrap_data, root, bond_num, -ClusterTree[root]);
        check_cross_wrap(wrap.X_wrap_roots, wrap.XY_wrap_roots, wrap.XY_wrap_data, root, bond_num, -ClusterTree[root]);
    }

    if( delta_x>1 && delta_y>1)
    {
        new_observed_wrap(wrap.XY_wrap_roots, wrap.XY_wrap_data, root, bond_num, -ClusterTree[root]);
    }

    if ( remove_from_chi && (wrap.XY_wrap_roots.empty() || wrap.XY_wrap_roots.find(root) == wrap.XY_wrap_roots.end()) ) remove_from_chi = 0;

    if (remove_from_chi)
    {
        const unsigned long long size = -ClusterTree[root];
        chi_num -= size*size;
        chi_den -= size;
    }

    check_master_equation(ClusterTree, wrap, chi_den, tot);
}


void check_master_equation(const std::vector<int>& ClusterTree, const WrapState& wrap, const unsigned long long& chi_den, const unsigned long long tot)
{
    unsigned long long PINFtot = 0;
    for (int root : wrap.XY_wrap_roots) PINFtot -= ClusterTree[root];
    if (PINFtot + chi_den != tot)
    {
        std::cout << "\nMaster equation is not verified!\n";
        std::cout << "Pinfs is given by: ";
        for (int root : wrap.XY_wrap_roots)  std::cout << -ClusterTree[root] << " ";
        std::cout << "\n";
        std::cout << "chi_den = " << chi_den << "\n";
        std::cout << "Their sum is " << PINFtot + chi_den << " given tot = " << tot << std::endl;
        exit(1);
    }
}


void measurements(const std::vector<int>& NZ, const std::vector<int>& RNZ, OrderParam& CP_Pinf, OrderParam& RP_Pinf,
                  LightOrderParam& CP_chi, LightOrderParam& RP_chi, Scalars& scalars)
{
    long double CP_chi_val, RP_chi_val;

    if (scalars.CP_chi_den>0) { CP_chi_val = (long double) scalars.CP_chi_num / (long double) scalars.CP_chi_den; } else { CP_chi_val = 0.0; }
    if (CP_chi_val > scalars.CP_chi_max) scalars.CP_chi_max = CP_chi_val;

    if (scalars.RP_chi_den>0) { RP_chi_val = (long double) scalars.RP_chi_num / (long double) scalars.RP_chi_den; } else { RP_chi_val = 0.0; }
    if (RP_chi_val > scalars.RP_chi_max) scalars.RP_chi_max = RP_chi_val;

    if (CP_Pinf.start <= (size_t)scalars.m && (size_t)scalars.m <= CP_Pinf.stop)
    {
        if (scalars.CPwrap.XY_wrap_roots.empty()) { update_OP(CP_Pinf, 0.0, 1.0, scalars.m - (int)CP_Pinf.start); }
        else
        {
            long double CP_Pinf_val = 0.0;
            for (size_t root : scalars.CPwrap.XY_wrap_roots) CP_Pinf_val += (long double)(-NZ[root]);
            update_OP(CP_Pinf, CP_Pinf_val, (long double)scalars.N, scalars.m - (int)CP_Pinf.start);
        }

        if (scalars.CP_chi_den>0) { update_OP_light(CP_chi, (long double)scalars.CP_chi_num, (long double)scalars.CP_chi_den, scalars.m - (int)CP_chi.start); }
        else { update_OP_light(CP_chi, 0.0, 1.0, scalars.m - (int)CP_chi.start); }
    }

    if (RP_Pinf.start <= (size_t)scalars.m && (size_t)scalars.m <= RP_Pinf.stop)
    {
        if (scalars.RPwrap.XY_wrap_roots.empty()) { update_OP(RP_Pinf, 0.0, 1.0, scalars.m - (int)RP_Pinf.start); }
        else
        {
            long double RP_Pinf_val = 0.0;
            for (size_t root : scalars.RPwrap.XY_wrap_roots) RP_Pinf_val += (long double)(-RNZ[root]);
            update_OP(RP_Pinf, RP_Pinf_val, (long double)scalars.M, scalars.m - (int)RP_Pinf.start);
        }

        if (scalars.RP_chi_den>0) { update_OP_light(RP_chi, (long double)scalars.RP_chi_num, (long double)scalars.RP_chi_den, scalars.m - (int)RP_chi.start); }
        else { update_OP_light(RP_chi, 0.0, 1.0, scalars.m - (int)RP_chi.start); }
    }
}


void find_relevant_wraps(const std::vector<WrappingTriplet>& x_wrap, const std::vector<WrappingTriplet>& xy_wrap, long int& x_first_bond, long int& xy_first_bond)
{
  x_first_bond = 1e12;
  xy_first_bond = 1e12;

  for (size_t i=0; i<x_wrap.size(); ++i)
  {
      const long int x_bond = x_wrap[i].bond_num;
      if (x_bond < x_first_bond) x_first_bond = x_bond;
  }

  for (size_t i=0; i<xy_wrap.size(); ++i)
  {
      const long int xy_bond = xy_wrap[i].bond_num;
      if (xy_bond < xy_first_bond) xy_first_bond = xy_bond;
  }

  if (x_first_bond == 1e12 || x_first_bond > xy_first_bond) x_first_bond = xy_first_bond;
}
