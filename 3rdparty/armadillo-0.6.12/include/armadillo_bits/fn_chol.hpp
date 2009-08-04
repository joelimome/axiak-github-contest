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


//! \addtogroup fn_chol
//! @{


template<typename eT, typename T1>
inline
bool
chol(Mat<eT>& out, const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(X.get_ref());
  arma_debug_check( !tmp.M.is_square(), "chol(): given matrix is not square");
  
  return auxlib::chol(out, tmp.M);
  }



template<typename eT, typename T1>
inline
Mat<eT>
chol(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  
  const bool ok = chol(out, X);
  if(ok == false)
    {
    arma_print("chol(): matrix factorisation failed");
    }
  
  return out;
  }


//! @}
