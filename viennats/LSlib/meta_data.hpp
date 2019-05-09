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
  class MetaData{

  public:
    std::vector< vec<double,D> > alpha;
    std::deque<bool> calculated; //NOTE vector<bool> is not used here, because it is not a STL container

    unsigned output_cnt;

    //TODO More elegant solution; maybe use generic static function
    void reinit_alpha(int n){
      alpha.resize(n);

      //fill all vectors with 0 entries
      for(unsigned alpha_idx = 0; alpha_idx < alpha.size(); ++alpha_idx){
        for(int i = 0; i < D; ++i)
          alpha[alpha_idx][i] = 0.0;
      }
    }

    void reinit_calculated(int n){
      calculated.resize(n);

      for(unsigned i = 0; i < calculated.size(); ++i){
          calculated[i] = false;
      }
    }

    MetaData():output_cnt(0){}

  };
}


#endif
