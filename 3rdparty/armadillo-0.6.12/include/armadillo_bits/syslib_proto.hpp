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


//! \addtogroup syslib
//! @{


class syslib
  {
  public:
  
  template<typename eT>
  arma_inline
  static
  void
  copy_elem(eT* dest, const eT* src, const u32 n_elem)
    {
    if( n_elem <= (128/sizeof(eT)) )
      {
      for(u32 i=0; i<n_elem; ++i)
        {
        dest[i] = src[i];
        }
      }
    else
      {
      std::memcpy(dest, src, n_elem*sizeof(eT));
      }
  
    }
  
  };


//! @}
