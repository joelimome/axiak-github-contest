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


//! \addtogroup fn_det
//! @{

template<typename T1> inline typename T1::elem_type det(const Base<typename T1::elem_type,T1>& X);
template<typename T1> inline typename T1::elem_type det(const Op<T1, op_diagmat>& X);

template<typename eT> inline eT det(const Op<Mat<eT>, op_diagmat_vec>& X);

template<typename T1, typename T2> inline typename T1::elem_type det(const Glue<T1, T2, glue_times>& X);

template<typename T1> inline typename T1::elem_type det(const Op<T1,op_inv>& in);
template<typename T1> inline typename T1::elem_type det(const Op<T1,op_trans>& in);



//! determinant of mat
template<typename T1>
inline
typename T1::elem_type
det(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> A_tmp(X.get_ref());
  const Mat<eT>& A = A_tmp.M;
  
  arma_debug_check( !A.is_square(), "det(): matrix must be square" );
  
  return auxlib::det(A);
  }



//! determinant of diagmat(mat)
template<typename T1>
inline
typename T1::elem_type
det(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();

  const unwrap<T1> A_tmp(X.m);

  typedef typename T1::elem_type eT;
  const Mat<eT>& A = A_tmp.M;

  arma_debug_check( (A.n_elem == 0), "det(): empty matrix");
  arma_debug_check( !A.is_square(), "det(): incompatible dimensions for diagmat operation" );

  eT val = A.at(0,0);
  
  for(u32 i=1; i<A.n_rows; ++i)
    {
    val *= A.at(i,i);
    }
  
  return val;
  }



//! determinant of diagmat(colvec or rowvec)
template<typename eT>
inline
eT
det(const Op<Mat<eT>, op_diagmat_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& A = X.m;
  
  arma_debug_check( (A.n_elem == 0), "det(): empty matrix");
  arma_debug_check( !A.is_vec(), "det_diagvec(): internal error: can't interpret as a vector" );

  eT val = A.mem[0];
  
  for(u32 i=1; i<A.n_elem; ++i)
    {
    val *= A.mem[i];
    }
  
  return val;
  }



//! determinant of A*B, avoiding the times operation if A and B are square matrices with the same dimensions
template<typename T1, typename T2>
inline
typename T1::elem_type
det(const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  if( (A.n_rows == A.n_cols) && (A.n_rows == B.n_rows) && (A.n_cols == B.n_cols) )
    {
    return det(A) * det(B);
    }
  else
    {
    return det(Mat<eT>(X));
    }
  
  }



//! determinant of inv(A), without doing the inverse operation
template<typename T1>
inline
typename T1::elem_type
det(const Op<T1,op_inv>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  isnt_fltpt<eT>::check();
  
  eT tmp = det(in.m);
  arma_debug_warn( (tmp == eT(0)), "det(): warning: determinant is zero" );
  
  return eT(1) / tmp;
  }



//! determinant of trans(A)
template<typename T1>
inline
typename T1::elem_type
det(const Op<T1,op_trans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;

  return det(X);
  }



//! @}
