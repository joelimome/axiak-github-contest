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


//! \addtogroup fn_trace
//! @{


//! Immediate trace (sum of diagonal elements) of a square dense matrix
template<typename T1>
inline
typename T1::elem_type
trace(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> A_tmp(X.get_ref());
  const Mat<eT>& A = A_tmp.M;

  arma_debug_check( !A.is_square(), "trace(): matrix must be square" );
  
  eT val = eT(0);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    val += A.at(i,i);
    }
  
  return val;
  }



//! \brief
//! Immediate trace (sum of diagonal elements) of A + B.
//! A and B must be square and have the same dimensions.
template<typename T1, typename T2>
inline
typename T1::elem_type
trace(const Glue<T1,T2,glue_plus>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols), "trace(): incompatible matrix dimensions");
  arma_debug_check( !A.is_square(), "trace(): matrices must be square");
  
  return trace(A) + trace(B);
  }



//! \brief
//! Immediate trace (sum of diagonal elements) of A - B.
//! A and B must be square and have the same dimensions.
template<typename T1, typename T2>
inline
typename T1::elem_type
trace(const Glue<T1,T2,glue_minus>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols), "trace(): incompatible matrix dimensions");
  arma_debug_check( !A.is_square(), "trace(): matrices must be square");
  
  return trace(A) - trace(B);
  }



//! \brief
//! Immediate trace (sum of diagonal elements) of A % B (where % is the element-wise multiplication operator).
//! A and B must be square and have the same dimensions.
template<typename T1, typename T2>
inline
typename T1::elem_type
trace(const Glue<T1,T2,glue_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();

  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols), "trace(): incompatible matrix dimensions" );
  arma_debug_check( !A.is_square(), "trace(): matrices must be square" );
  
  eT val = eT(0);
  for(u32 i=0; i<A.n_rows; ++i)
    {
    val += A.at(i,i) * B.at(i,i);
    }
  
  return val;
  }



//! \brief
//! trace (sum of diagonal elements) of k * T1,
//! where k is a scalar and T1 is converted to a dense matrix.
template<typename T1>
inline
typename T1::elem_type
trace(const Op<T1,op_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  return trace(X) * in.aux;
  }



//! trace (sum of diagonal elements) of a diagonal matrix
template<typename eT>
inline
eT
trace(const Op<Mat<eT>, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  return trace(X.m);
  }



template<typename eT>
inline
eT
trace(const Op<Mat<eT>, op_diagmat_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& A = X.m;
  arma_debug_check( !A.is_vec(), "trace(): internal error: can't interpret as a vector" );
  
  
  return accu(X.m);
  }


//! @}
