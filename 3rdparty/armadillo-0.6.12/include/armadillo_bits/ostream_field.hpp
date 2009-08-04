// Copyright (C) 2009 NICTA
// 
// Authors:
// - Conrad Sanderson (conradsand at ieee dot org)
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup ostream
//! @{


//! Print the contents of a field to the specified stream
//! Assumes type T1 can be printed, i.e. T1 has std::ostream& operator<< (std::ostream&, const T1&) 

template<typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const field<T1>& X)
  {
  arma_extra_debug_sigprint();
  
  for(u32 col=0; col<X.n_cols; ++col)
    {
    o << "[field column " << col << ']' << '\n'; 
    for(u32 row=0; row<X.n_rows; ++row)
      {
      o << X.at(row,col) << '\n';
      }
    
    o << '\n';
    }
  
  o.flush();
  
  return o;
  }



//! Print the contents of a subfield to the specified stream
//! Assumes type T1 can be printed, i.e. T1 has std::ostream& operator<< (std::ostream&, const T1&) 

template<typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const subview_field<T1>& X)
  {
  arma_extra_debug_sigprint();
  
  for(u32 col=0; col<X.n_cols; ++col)
    {
    for(u32 row=0; row<X.n_rows; ++row)
      {
      o << X.at(row,col) << '\n';
      }
    
    o << '\n';
    }
  
  o.flush();
  
  return o;
  }



//! @}
