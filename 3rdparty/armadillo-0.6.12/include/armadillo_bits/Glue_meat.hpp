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


//! \addtogroup Glue
//! @{


template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::Glue(const T1& in_A, const T2& in_B)
  : A(in_A)
  , B(in_B)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  }



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::~Glue()
  {
  arma_extra_debug_sigprint();
  }


//! @}
