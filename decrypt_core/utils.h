// Copyright 2020 Arnaud Becheler    <arnaud.becheler@gmail.com>

/***********************************************************************                                                                         *
* This program is free software; you can redistribute it and/or modify *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation; either version 2 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
************************************************************************/

#ifndef __DECRYPT_UTILS_H_INCLUDED__
#define __DECRYPT_UTILS_H_INCLUDED__

#include "quetzal.h"

namespace decrypt::utils{

  // used in dispersal kernel definition
  template<typename T>
  auto make_neighboring_cells_functor(T const& landscape){

    using coord_type = quetzal::geography::GeographicCoordinates;
    return [landscape](coord_type const& x0){
      auto res = landscape.get().resolution();

      std::vector<coord_type> v;

      coord_type x1(x0);
      x1.lon() += res.lon();
      if(landscape.get().is_in_spatial_extent(x1)){
        v.push_back(landscape.get().reproject_to_centroid(x1));
      }

      coord_type x2(x0);
      x2.lat() += res.lat();
      if(landscape.get().is_in_spatial_extent(x2)){
        v.push_back(landscape.get().reproject_to_centroid(x2));
      }

      coord_type x3(x0);
      x3.lon() -= res.lon();
      if(landscape.get().is_in_spatial_extent(x3)){
        v.push_back(landscape.get().reproject_to_centroid(x3));
      }

      coord_type x4(x0);
      x4.lat() -= res.lat() ;
      if(landscape.get().is_in_spatial_extent(x4)){
        v.push_back(landscape.get().reproject_to_centroid(x4));
      }
      return v;
    };
  };

  // Class for representing sampled individual gene copies
  class GeneCopy
  {
  private:
    using coord_type = quetzal::geography::GeographicCoordinates;
    std::string m_id;
    coord_type m_x;
  public:
    GeneCopy(std::string const& id, coord_type const& x): m_id(id), m_x(x){}
    auto id() const {return m_id;}
    auto x() const {return m_x;}
  };
}

#endif
