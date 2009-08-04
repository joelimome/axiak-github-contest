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


//! \addtogroup fn_accu
//! @{



//! accumulate the elements of a matrix
template<typename T1>
inline
typename T1::elem_type
accu(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(X.get_ref());
  const Mat<eT>& A = tmp.M;
  
  const u32 A_n_elem = A.n_elem;
  const eT* A_mem    = A.mem;
  
  eT val = eT(0);
  
  for(u32 i=0; i<A_n_elem; ++i)
    {
    val += A_mem[i];
    }
  
  return val;
  }



//! sum of values along the main diagonal
template<typename T1>
inline
typename T1::elem_type
accu(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(X.m);
  const Mat<eT>& A = tmp.M;
  
  arma_debug_check( !A.is_square(), "accu(): sum of diagonal values of a non-square matrix requested" );
  
  eT acc = eT(0);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    acc += A.at(i,i);
    }
  
  return acc;
  }



template<typename eT>
inline
eT
accu(const Op<Mat<eT>, op_diagmat_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& A = X.m;
  arma_debug_check( !A.is_vec(), "accu(): internal error: expected a vector" );
  
  return accu(A);
  }



//! sum of squares
template<typename T1>
inline
typename T1::elem_type
accu(const Op<T1, op_square>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  const u32 A_n_elem = A.n_elem;
  const eT* A_mem    = A.mem;
  
  eT acc = eT(0);
  
  for(u32 i=0; i<A_n_elem; ++i)
    {
    const eT val = A_mem[i];
    acc += val*val;
    }
  
  return acc;
  }



//! sum of square roots
template<typename T1>
inline
typename T1::elem_type
accu(const Op<T1, op_sqrt>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  const u32 A_n_elem = A.n_elem;
  const eT* A_mem    = A.mem;
  
  eT acc = eT(0);
  for(u32 i=0; i<A_n_elem; ++i)
    {
    acc += std::sqrt(A_mem[i]);
    }
  
  return acc;
  }



//! sum of squares of differences
template<typename T1, typename T2>
inline
typename T1::elem_type
accu(const Op< Glue<T1,T2, glue_minus>, op_square>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A,B, "accu()");
  
  const u32 n_elem = A.n_elem;
  
  eT acc = eT(0);
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT val = A.mem[i] - B.mem[i];
    acc += val*val;
    }
  
  return acc;
  }



//! accumulate the elements of a subview (submatrix)
template<typename eT>
inline
eT
accu(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  eT val = eT(0);
  for(u32 col=0; col<X.n_cols; ++col)
    {
    const eT* coldata = X.colptr(col);
    
    for(u32 row=0; row<X.n_rows; ++row)
      {
      val += coldata[row];
      }
    
    }
  
  return val;
  }



//! accumulate the elements of a diagview
template<typename eT>
inline
eT
accu(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();  
  
  const u32 n_elem = X.n_elem;
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += X[i];
    }
  
  return val;
  }



//! accumulate the result of A % B, where % is the Schur product (element-wise multiplication)
template<typename eT>
inline
eT
accu_schur(const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A,B, "accu()");
 
  const eT* const A_mem = A.mem;
  const eT* const B_mem = B.mem;
  
  const u32 n_elem = A.n_elem;
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += A_mem[i] * B_mem[i];
    }
  
  return val;
  }



//! accumulate the result of A % B, where % is the Schur product (element-wise multiplication)
template<typename eT>
inline
eT
accu(const Glue<Mat<eT>,Mat<eT>,glue_schur>& X)
  {
  return accu_schur(X.A, X.B);
  }



//! accumulate the result of A % B % C, where % is the Schur product (element-wise multiplication)
template<typename eT>
inline
eT
accu(const Glue<Glue<Mat<eT>,Mat<eT>,glue_schur>,Mat<eT>,glue_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& A = X.A.A;
  const Mat<eT>& B = X.A.B;
  const Mat<eT>& C = X.B;
  
  arma_debug_assert_same_size(A,B, "accu()");
  arma_debug_assert_same_size(A,C, "accu()");
  
  const eT* const A_mem = A.mem;
  const eT* const B_mem = B.mem;
  const eT* const C_mem = C.mem;
  
  const u32 n_elem = A.n_elem;
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += A_mem[i] * B_mem[i] * C_mem[i];
    }
    
  return val;
  }



//! \brief
//! accumulate the result of T1 % T2 
//! where % is the Schur product (element-wise multiplication),
//! while T1 and T2 can be 'mat', 'rowvec', 'colvec', 'Op', 'Glue'
    
template<typename T1, typename T2>
inline
typename T1::elem_type
accu(const Glue<T1,T2,glue_schur>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();

  typedef typename T1::elem_type eT;

  const u32 N_mat = 1 + depth_lhs< glue_schur, Glue<T1,T2,glue_schur> >::num;
  arma_extra_debug_print( arma_boost::format("N_mat = %d") % N_mat );

  if(N_mat == 2)
    {
    return accu_schur(Mat<eT>(X.A), Mat<eT>(X.B));
    }
  else
    {
    const Mat<eT>* ptrs[N_mat];
    bool            del[N_mat];
  
    mat_ptrs<glue_schur, Glue<T1,T2,glue_schur> >::get_ptrs(ptrs, del, X);
  
    for(u32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format("ptrs[%d] = %x") % i % ptrs[i] );
    for(u32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format(" del[%d] = %d") % i %  del[i] );
  
    const Mat<eT>& tmp_mat = *(ptrs[0]);
    
    for(u32 i=1; i<N_mat; ++i)
      {
      arma_debug_assert_same_size(tmp_mat, *(ptrs[i]), "accu()");
      }
    
    // const u32 n_rows = ptrs[0]->n_rows;
    // const u32 n_cols = ptrs[0]->n_cols;
    
    eT val = eT(0);
    
    const u32 n_elem = ptrs[0]->n_elem;
    
    for(u32 j=0; j<n_elem; ++j)
      {
      eT tmp = ptrs[0]->mem[j];
    
      for(u32 i=1; i<N_mat; ++i)
        {
        tmp *= ptrs[i]->mem[j];
        }
    
      val += tmp;
      }
    
    
    for(u32 i=0; i<N_mat; ++i)
      {
      if(del[i] == true)
        {
        arma_extra_debug_print( arma_boost::format("delete mat_ptr[%d]") % i );
        delete ptrs[i];
        }
      }
    
    return val;    
    }
  }


//! accumulate the result of submatrix % matrix, where % is the Schur product (element-wise multiplication)
template<typename eT>
inline
eT
accu(const Glue<subview<eT>,Mat<eT>,glue_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(X.A, X.B, "accu()");
  
  const Mat<eT>& A = X.A.m;
  const Mat<eT>& B = X.B;
  
  const u32 A_sub_n_rows = X.A.n_rows;
  const u32 A_sub_n_cols = X.A.n_cols;
  
  const u32 A_aux_row1 = X.A.aux_row1;
  const u32 A_aux_col1 = X.A.aux_col1;
  
  
  eT val = eT(0);
    
  for(u32 col = 0; col<A_sub_n_cols; ++col)
    {
    const u32 col_mod = A_aux_col1 + col;
    
    for(u32 row = 0; row<A_sub_n_rows; ++row)
      {
      const u32 row_mod = A_aux_row1 + row;
      
      val += A.at(row_mod, col_mod) * B.at(row,col);
      }
    
    }
  
  return val;
  }



//! accumulate the result of matrix % submatrix, where % is the Schur product (element-wise multiplication)
template<typename eT>
inline
eT
accu(const Glue<Mat<eT>,subview<eT>,glue_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(X.A, X.B, "accu()");
  
  const Mat<eT>& A = X.A;
  const Mat<eT>& B = X.B.m;
  
  // const u32 B_sub_n_rows = X.B.n_rows;
  // const u32 B_sub_n_cols = X.B.n_cols;
  
  const u32 B_aux_row1 = X.B.aux_row1;
  const u32 B_aux_col1 = X.B.aux_col1;
  
  
  eT val = eT(0);
    
  for(u32 col = 0; col<A.n_cols; ++col)
    {
    const u32 col_mod = B_aux_col1 + col;
    
    for(u32 row = 0; row<A.n_rows; ++row)
      {
      const u32 row_mod = B_aux_row1 + row;
      
      val += A.at(row, col) * B.at(row_mod, col_mod);
      }
    
    }
  
  return val;
  }



//! accumulate the result of submatrix % submatrix, where % is the Schur product (element-wise multiplication)
template<typename eT>
inline
eT
accu(const Glue<subview<eT>,subview<eT>,glue_schur>& X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(X.A, X.B, "accu()");
  
  const Mat<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m;
  
  const u32 A_sub_n_rows = X.A.n_rows;
  const u32 A_sub_n_cols = X.A.n_cols;
  
  // const u32 B_sub_n_rows = X.B.n_rows;
  // const u32 B_sub_n_cols = X.B.n_cols;
  
  const u32 A_aux_row1 = X.A.aux_row1;
  const u32 A_aux_col1 = X.A.aux_col1;
  
  const u32 B_aux_row1 = X.B.aux_row1;
  const u32 B_aux_col1 = X.B.aux_col1;
  
  
  eT val = eT(0);
    
  for(u32 col = 0; col<A_sub_n_cols; ++col)
    {
    const u32 A_col_mod = A_aux_col1 + col;
    const u32 B_col_mod = B_aux_col1 + col;
    
    for(u32 row = 0; row<A_sub_n_rows; ++row)
      {
      const u32 A_row_mod = A_aux_row1 + row;
      const u32 B_row_mod = B_aux_row1 + row;
      
      val += A.at(A_row_mod, A_col_mod) * B.at(B_row_mod, B_col_mod);
      }
    
    }
  
  return val;
  }



//! @}
