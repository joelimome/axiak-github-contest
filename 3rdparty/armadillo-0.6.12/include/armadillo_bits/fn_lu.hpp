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


//! \addtogroup fn_lu
//! @{


//! immediate lower upper decomposition
template<typename eT, typename T1>
inline
void
lu(Mat<eT>& L, Mat<eT>& U, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(X.get_ref());

  auxlib::lu(L, U, tmp.M);
  }



//! immediate lower upper decomposition, also providing the permutation matrix
template<typename eT, typename T1>
inline
void
lu(Mat<eT>& L, Mat<eT>& U, Mat<eT>& P, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(X.get_ref());
  
  auxlib::lu(L, U, P, tmp.M);
  }


//! @}
