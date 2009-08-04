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


//! \addtogroup op_inv
//! @{


//! immediate inverse of a matrix, storing the result in a dense matrix
template<typename eT>
inline
void
op_inv::apply(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  // no need to check for aliasing, due to:
  // - auxlib::inv() copies A to out before inversion
  // - for 2x2 and 3x3 matrices the code is alias safe
  
  arma_debug_check( !A.is_square(), "op_inv::apply(): matrix must be square" );
  
  if(&out != &A)
    {
    auxlib::inv_noalias(out, A);
    }
  else
    {
    auxlib::inv_inplace(out);
    }
  
  }



//! immediate inverse of T1, storing the result in a dense matrix
template<typename T1>
inline
void
op_inv::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_inv>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.m);
  
  typedef typename T1::elem_type eT;
  const Mat<eT>& X = tmp.M;
  
  op_inv::apply(out, X);
  }



//! inverse of diagmat(mat)
template<typename T1>
inline
void
op_inv::apply(Mat<typename T1::elem_type>& out, const Op< Op<T1,op_diagmat>, op_inv>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> X_tmp(in.m.m);
  
  typedef typename T1::elem_type eT;
  const Mat<eT>& X = X_tmp.M;
  
  arma_debug_check( !X.is_square(), "op_inv::apply(): matrix must be square" );
  
  if(&out != &X)
    {
    out.zeros(X.n_rows, X.n_rows);
    
    for(u32 i=0; i<X.n_rows; ++i)
      {
      out.at(i,i) = 1.0 / X.at(i,i);
      }
    }
  else
    {
    podarray<eT> tmp(X.n_rows);
    
    for(u32 i=0; i<X.n_rows; ++i)
      {
      tmp[i] = X.at(i,i);
      }
      
    out.zeros(X.n_rows, X.n_rows);
    
    for(u32 i=0; i<X.n_rows; ++i)
      {
      out.at(i,i) = eT(1) / tmp.mem[i];
      }
    
    }
  
  }



template<typename eT>
inline
void
op_inv::apply_diagvec(Mat<eT>& out, const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();

  arma_debug_check( !X.is_vec(), "op_inv::apply_diagvec(): internal error: can't interpret as a vector");
  
  if(&out != &X)
    {
    out.zeros(X.n_elem, X.n_elem);
    
    for(u32 i=0; i<X.n_elem; ++i)
      {
      out.at(i,i) = eT(1) / X.mem[i];
      }
    }
  else
    {
    podarray<eT> tmp(X.n_elem);
    
    for(u32 i=0; i<X.n_elem; ++i)
      {
      tmp[i] = X.mem[i];
      }
      
    out.zeros(X.n_elem, X.n_elem);
    
    for(u32 i=0; i<X.n_elem; ++i)
      {
      out.at(i,i) = eT(1) / tmp.mem[i];
      }
    
    }

  }



//! inverse of diagmat(colvec or rowvec)
template<typename eT>
inline
void
op_inv::apply(Mat<eT>& out, const Op< Op<Mat<eT>,op_diagmat_vec>, op_inv>& in)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& X = in.m.m;
  op_inv::apply_diagvec(out, X);
  }



//! @}
