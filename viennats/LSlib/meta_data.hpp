#ifndef META_DATA_HPP_
#define META_DATA_HPP_

/* =========================================================================
   Copyright (c)    2008-2019, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


#include <vector>
#include <deque>
#include "levelset.hpp"

namespace lvlset {


  template<int D>
  class ExpandedPointsData{
    public:
      std::vector< vec<double,D> > vec_values;
      std::deque<bool> flags; //NOTE vector<bool> is not used here, because it is not a STL container

      void reinit(int n){
        vec_values.resize(n);

        //fill all vectors with 0 entries
        for(unsigned values_idx = 0; values_idx < vec_values.size(); ++values_idx){
          for(int i = 0; i < D; ++i)
            vec_values[values_idx][i] = 0.0;
        }

        flags.resize(n);

        for(unsigned i = 0; i < flags.size(); ++i){
            flags[i] = false;
        }
      }
  };

  template<int D>
  class ActivePointsData{
    public:
      std::vector<double > scalar_values;

      void reinit(int n){
        scalar_values.resize(n);

        for(unsigned i = 0; i < scalar_values.size(); ++i){
            scalar_values[i] = 0.0;
        }
      }
  };

  template<int D>
  class MetaData{
    public:
      ExpandedPointsData<D> expanded_pt_data;
      ActivePointsData<D> active_pt_data;
      unsigned output_cnt;


      MetaData():output_cnt(0){}

  };
}


#endif
