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


//! \addtogroup op_sum
//! @{

//! \brief
//! Immediate sum of elements of a matrix along a specified dimension (either rows or columns).
//! The result is stored in a dense matrix that has either one column or one row.
//! See the sum() function for more details.
template<typename T1>
inline
void
op_sum::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sum>& in)
  {
  arma_extra_debug_sigprint();
  
  const u32 dim = in.aux_u32_a;
  arma_debug_check( (dim > 1), "op_sum::apply(): incorrect usage. dim must be 0 or 1");
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp(in.m, out);
  const Mat<eT>& X = tmp.M;
  
  if(dim == 0)  // column-wise sum
    {
    out.set_size(1, X.n_cols);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      eT tmp_val = eT(0);
      for(u32 row=0; row<X.n_rows; ++row)
        {
        tmp_val += X.at(row,col);
        }
    
      out.at(0,col) = tmp_val;
      }
    }
  else
  if(dim == 1)  // row-wise sum
    {
    out.set_size(X.n_rows, 1);
    
    for(u32 row=0; row<X.n_rows; ++row)
      {
      eT tmp_val = eT(0);
      for(u32 col=0; col<X.n_cols; ++col)
        {
        tmp_val += X.at(row,col);
        }
    
      out.at(row,0) = tmp_val;
      }
    
    }
  
  }

//! @}
