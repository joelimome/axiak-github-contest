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



//! \addtogroup op_neg
//! @{

//! Negate all element of a matrix and store the result in a dense matrix
template<typename T1>
inline
void
op_neg::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_neg> &in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // no alias problems
  out.set_size(X.n_rows, X.n_cols);
  
  const eT* X_mem = X.mem;
  eT* out_mem = out.memptr();
  
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out_mem[i] = -X_mem[i];
    }
    
  }

//! @}
