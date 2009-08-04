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


//! \addtogroup op_diagmat
//! @{



template<typename eT>
inline
void
op_diagmat::zero_offdiag(Mat<eT>& X)
  {
  for(u32 col=0; col<X.n_cols; ++col)
    {
    eT* colmem = X.colptr(col);
    
    // above the diagonal
    
    for(u32 row=0; row<col; ++row)
      {
      colmem[row] = eT(0);
      }
    
    // below the diagonal
    
    for(u32 row=col+1; row<X.n_rows; ++row)
      {
      colmem[row] = eT(0);
      }
    }
  }



template<typename T1>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.m);
  
  typedef typename T1::elem_type eT;
  const Mat<eT>& X = tmp.M;
  
  arma_debug_check( !X.is_square(), "diagmat(): matrix must be square" );
  
  if(&out != &X)
    {
    out.zeros(X.n_rows, X.n_rows);
    
    for(u32 i=0; i<X.n_rows; ++i)
      {
      out.at(i,i) = X.at(i,i);
      }
    }
  else
    {
    op_diagmat::zero_offdiag(out);
    }
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_div>, op_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (!A.is_square() || !B.is_square()), "diagmat(): matrices must be square" );
  arma_debug_assert_same_size(A, B, "matrix element-wise division");
    
  // not using out.zeros() as 'out' might be an alias of A and/or B.
  // the off-diagonal elements are zeroed at the end
  out.set_size(A.n_rows, A.n_rows);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = A.at(i,i) / B.at(i,i);
    }
  
  op_diagmat::zero_offdiag(out);
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_minus>, op_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (!A.is_square() || !B.is_square()), "diagmat(): matrices must be square" );
  arma_debug_assert_same_size(A, B, "matrix subtraction");
    
  // not using out.zeros() as 'out' might be an alias of A and/or B.
  // the off-diagonal elements are zeroed at the end
  out.set_size(A.n_rows, A.n_rows);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = A.at(i,i) - B.at(i,i);
    }
  
  op_diagmat::zero_offdiag(out);
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_plus>, op_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (!A.is_square() || !B.is_square()), "diagmat(): matrices must be square" );
  arma_debug_assert_same_size(A, B, "matrix addition");
    
  // not using out.zeros() as 'out' might be an alias of A and/or B.
  // the off-diagonal elements are zeroed at the end
  out.set_size(A.n_rows, A.n_rows);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = A.at(i,i) + B.at(i,i);
    }
  
  op_diagmat::zero_offdiag(out);
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_schur>, op_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (!A.is_square() || !B.is_square()), "diagmat(): matrices must be square" );
  arma_debug_assert_same_size(A, B, "matrix schur product");
    
  // not using out.zeros() as 'out' might be an alias of A and/or B.
  // the off-diagonal elements are zeroed at the end
  out.set_size(A.n_rows, A.n_rows);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = A.at(i,i) * B.at(i,i);
    }
  
  op_diagmat::zero_offdiag(out);
  }



template<typename T1, typename T2>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_times>, op_diagmat>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap_check<T1> tmp1(in.m.A, out);
  const unwrap_check<T2> tmp2(in.m.B, out);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.n_rows != B.n_cols), "diagmat(): result of multiplication is not square" );
  arma_debug_assert_mul_size(A, B, "matrix multiplication");
    
  // out is cleared here, as we've made sure that A and B are not aliases of 'out'
  out.zeros(A.n_rows, A.n_rows);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    const eT* B_colmem = B.colptr(i);
    
    eT val = eT(0);
    for(u32 j=0; j<A.n_cols; ++j)
      {
      val += A.at(i,j) * B_colmem[j];
      }
    
    out.at(i,i) = val;
    }
  }



//
//
//


template<typename T1>
inline
void
op_diagmat_vec::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagmat_vec>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp(in.m);
  
  typedef typename T1::elem_type eT;
  const Mat<eT>& X = tmp.M;
  
  if(&out != &X)
    {
    out.zeros(X.n_elem, X.n_elem);
    
    const eT* X_mem = X.mem;
    
    for(u32 i=0; i<X.n_elem; ++i)
      {
      out.at(i,i) = X_mem[i];
      }
    }
  else
    {
    podarray<eT> tmp_array(X.n_elem);
    
    for(u32 i=0; i<X.n_elem; ++i)
      {
      tmp_array[i] = X[i];
      }
    
    out.zeros(tmp_array.n_elem, tmp_array.n_elem);
    
    for(u32 i=0; i<tmp_array.n_elem; ++i)
      {
      out.at(i,i) = tmp_array[i];
      }
    
    }
  }



//! @}
