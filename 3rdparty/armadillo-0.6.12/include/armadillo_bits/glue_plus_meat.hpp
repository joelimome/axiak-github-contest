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


//! \addtogroup glue_plus
//! @{



template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  // no aliasing problem
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
    
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A_mem[i] + B_mem[i];
    }
    
  }



template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  arma_debug_assert_same_size(A, C, "matrix addition");
  
  // no aliasing problem
  out.set_size(A.n_rows, A.n_cols);
    
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  const eT* C_mem   = C.mem;
  
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = A_mem[i] + B_mem[i] + C_mem[i];
    }
    
  }



template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Glue<Mat<eT>,Mat<eT>,glue_plus>& X)
  {
  glue_plus::apply(out, X.A, X.B);
  }



template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Glue< Glue<Mat<eT>,Mat<eT>,glue_plus>, Mat<eT>, glue_plus>& X)
  {
  glue_plus::apply(out, X.A.A, X.A.B, X.B);
  }



template<typename T1, typename T2>
inline
void
glue_plus::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const u32 N_mat = 1 + depth_lhs< glue_plus, Glue<T1,T2,glue_plus> >::num;
  arma_extra_debug_print( arma_boost::format("N_mat = %d") % N_mat );

  if(N_mat == 2)
    {
    if(is_glue_times<T1>::value == true)
      {
      out = X.B;
      glue_plus::apply_inplace(out, X.A);
      }
    else
    if(is_glue_times<T2>::value == true)
      {
      out = X.A;
      glue_plus::apply_inplace(out, X.B);
      }
    else
      {
      const unwrap<T1> tmp1(X.A);
      const unwrap<T2> tmp2(X.B);
      
      glue_plus::apply(out, tmp1.M, tmp2.M);
      }
    }
  else
    {
    const Mat<eT>* ptrs[N_mat];
    bool            del[N_mat];

    mat_ptrs<glue_plus, Glue<T1,T2,glue_plus> >::get_ptrs(ptrs, del, X);

    for(u32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format("ptrs[%d] = %x") % i % ptrs[i] );
    for(u32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format(" del[%d] = %d") % i %  del[i] );

    const u32 n_rows = ptrs[0]->n_rows;
    const u32 n_cols = ptrs[0]->n_cols;
  
    const Mat<eT>& tmp_mat = *(ptrs[0]);
    
    for(u32 i=1; i<N_mat; ++i)
      {
      arma_debug_assert_same_size(tmp_mat, *(ptrs[i]), "matrix addition");
      }
  
  
    // no aliasing problem
    out.set_size(n_rows, n_cols);
    
    const u32 n_elem = ptrs[0]->n_elem;
    
    for(u32 j=0; j<n_elem; ++j)
      {
      eT acc = ptrs[0]->mem[j];
    
      for(u32 i=1; i < N_mat; ++i)
        {
        acc += ptrs[i]->mem[j];
        }
    
      out[j] = acc;
      }
    
    
    for(u32 i=0; i<N_mat; ++i)
      {
      if(del[i] == true)
        {
        arma_extra_debug_print( arma_boost::format("delete mat_ptr[%d]") % i );
        delete ptrs[i];
        }
      }
    }
  }



// possible aliasing cases:
// Q = Q + Q.row(0)  -> no problem  (aliasing has no effect or incompatible matrix dimensions).
//                      however, the only time the above will work is when Q has the same dimensions as Q.row(0),
//                      meaning that doing this addition operation is pretty silly
// Q = Q + R.row(0)  -> no problem
// Q = R + Q.row(0)  -> output Q is set to size of R, which may destroy input Q
//
// strategy:
// if the matrix from the second argument is an alias of the output matrix,
// make a proper matrix out of the second argument

template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Glue<Mat<eT>, subview<eT>, glue_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& orig_A = X.A;
  const Mat<eT>& orig_B = X.B.m;
  
  if( &out != &orig_B )
    {
    //const u32 sub_B_n_rows = X.B.n_rows;
    //const u32 sub_B_n_cols = X.B.n_cols;
    
    arma_debug_assert_same_size(X.A, X.B, "matrix addition");
      
    
    out.set_size(orig_A.n_rows, orig_A.n_cols);
    
    for(u32 col = 0; col<orig_A.n_cols; ++col)
      {
      const u32 B_col_mod = X.B.aux_col1 + col;
      
      for(u32 row = 0; row<orig_A.n_rows; ++row)
        {
        const u32 B_row_mod = X.B.aux_row1 + row;
        
        out.at(row,col) =  orig_A.at(row, col) + orig_B.at(B_row_mod, B_col_mod);
        }
      
      }
    
    }
  else
    {
    const Mat<eT> processed_B(X.B);  // create a matrix out of subview
    glue_plus::apply(out, orig_A, processed_B);
    }
     
  }


// possible aliasing cases:
// Q = Q.row(0) + Q  -> no problem (aliasing has no effect or incompatible matrix dimensions)
// Q = Q.row(0) + R  -> problem (output Q is set to size of Q.row(0) which may destroy input Q)
// Q = R.row(0) + Q  -> problem (output Q is set to size of R.row(0) which may destroy input Q)

template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Glue< subview<eT>, Mat<eT>, glue_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& orig_A = X.A.m;
  
  const unwrap_check< Mat<eT> > tmp(X.B, out);
  const Mat<eT>& orig_B = tmp.M;
  
  if( &out != &orig_A )
    {
    const u32 sub_A_n_rows = X.A.n_rows;
    const u32 sub_A_n_cols = X.A.n_cols;
    
    arma_debug_assert_same_size(X.A, X.B, "matrix addition");
      
    out.set_size(sub_A_n_rows, sub_A_n_cols);
    
    for(u32 col = 0; col<sub_A_n_cols; ++col)
      {
      const u32 A_col_mod = X.A.aux_col1 + col;
      
      for(u32 row = 0; row<sub_A_n_rows; ++row)
        {
        const u32 A_row_mod = X.A.aux_row1 + row;
        
        out.at(row,col) =  orig_A.at(A_row_mod, A_col_mod) + orig_B.at(row, col);
        }
      
      }
    }
  else
    {
    const Mat<eT> processed_A(X.A);
    glue_plus::apply(out, processed_A, orig_B);
    }
  
  
  }


// possible aliasing cases:
// Q = Q.row(0) + Q.row(0)  -> input Q is destroyed unless Q.row(0) has the same size as Q 
// Q = Q.row(0) + R.row(0)  -> input Q is destroyed unless Q.row(0) has the same size as Q 
// Q = R.row(0) + Q.row(0)  -> input Q is destroyed

template<typename eT>
inline
void
glue_plus::apply(Mat<eT>& out, const Glue< subview<eT>, subview<eT>, glue_plus>& X)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& orig_A = X.A.m;
  const Mat<eT>& orig_B = X.B.m;
  
  if( (&out != &orig_A) && (&out != &orig_B) )
    {
    const u32 sub_A_n_rows = X.A.n_rows;
    const u32 sub_A_n_cols = X.A.n_cols;
    
    //const u32 sub_B_n_rows = X.B.n_rows;
    //const u32 sub_B_n_cols = X.B.n_cols;
    
    arma_debug_assert_same_size(X.A, X.B, "matrix addition");
      
    out.set_size(sub_A_n_rows, sub_A_n_cols);
    
    for(u32 col = 0; col<sub_A_n_cols; ++col)
      {
      const u32 A_col_mod = X.A.aux_col1 + col;
      const u32 B_col_mod = X.B.aux_col1 + col;
      
      for(u32 row = 0; row<sub_A_n_rows; ++row)
        {
        const u32 A_row_mod = X.A.aux_row1 + row;
        const u32 B_row_mod = X.B.aux_row1 + row;
        
        out.at(row,col) =  orig_A.at(A_row_mod, A_col_mod) + orig_B.at(B_row_mod, B_col_mod);
        }
      
      }
    }
  else
    {
    const Mat<eT> processed_A(X.A);
    const Mat<eT> processed_B(X.B);
    
    glue_plus::apply(out, processed_A, processed_B);
    }
  }



template<typename eT>
inline void glue_plus::apply_inplace(Mat<eT>& out, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, B, "matrix addition");
  
  
        eT* out_mem = out.memptr();
  const eT* B_mem   = B.mem;
  
  const u32 n_elem  = B.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] += B_mem[i];
    }
  
  }



template<typename T1, typename op_type>
inline
void
glue_plus::apply_inplace(Mat<typename T1::elem_type>& out, const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT> tmp(X);
  glue_plus::apply(out, out, tmp);
  }
  


template<typename T1>
inline
void
glue_plus::apply_inplace(Mat<typename T1::elem_type>& out, const Op<T1, op_square>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(X.m);
  const Mat<eT>& B = tmp.M;
  
  arma_debug_assert_same_size(out, B, "matrix addition");
    
        eT* out_mem = out.memptr();
  const eT* B_mem   = B.mem;
  
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT tmp_val = B_mem[i];
    out_mem[i] += tmp_val*tmp_val;
    }
  }



template<typename T1, typename T2, typename glue_type>
inline
void
glue_plus::apply_inplace(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
    
  out = X + out;
  }



template<typename T1, typename T2>
inline
void
glue_plus::apply_inplace(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp1(X.A, out);
  const unwrap_check<T2> tmp2(X.B, out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "matrix multiplication");
  arma_debug_assert_same_size(out.n_rows, out.n_cols, A.n_rows, B.n_cols, "matrix addition");
  
  gemm<false,false,false,true>::apply(out, A, B, eT(1), eT(1));
  }



//
//
//



template<typename T1, typename T2>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue< Glue< T1, Col<typename T1::elem_type>, glue_times_vec>, T2, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.A);
  const unwrap<T2> tmp2(in.B);
  
  const Mat<eT>& A = tmp1.M;
  const Col<eT>& B = in.A.B;
  const Mat<eT>& C = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "matrix multiplication");
  arma_debug_assert_same_size(A.n_rows, B.n_cols, C.n_rows, C.n_cols, "matrix addition");
  
  if( (&out != &A) && (&out != &B) )
    {
    out = C;
    gemv<false,false,true>::apply(out.memptr(), A, B.mem, eT(1), eT(1));
    }
  else
    {
    const unwrap_check< Mat<eT> > tmpA(A,out);
    const unwrap_check< Col<eT> > tmpB(B,out);
    
    const Mat<eT>& A_safe = tmpA.M;
    const Col<eT>& B_safe = tmpB.M;
    
    out = C;
    gemv<false,false,true>::apply(out.memptr(), A_safe, B_safe.mem, eT(1), eT(1));
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue< Glue< Row<typename T1::elem_type>, T1, glue_times_vec>, T2, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.B);
  const unwrap<T2> tmp2(in.B);
  
  const Row<eT>& A = in.A.A;
  const Mat<eT>& B = tmp1.M;
  const Mat<eT>& C = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "matrix multiplication");
  arma_debug_assert_same_size(A.n_rows, B.n_cols, C.n_rows, C.n_cols, "matrix addition");
  
  if( (&out != &A) && (&out != &B) )
    {
    out = C;
    gemv<true,false,true>::apply(out.memptr(), B, A.mem, eT(1), eT(1));
    }
  else
    {
    const unwrap_check< Mat<eT> > tmpA(A,out);
    const unwrap_check< Mat<eT> > tmpB(B,out);
    
    const Mat<eT>& A_safe = tmpA.M;
    const Mat<eT>& B_safe = tmpB.M;
    
    out = C;
    gemv<true,false,true>::apply(out.memptr(), B_safe, A_safe.mem, eT(1), eT(1));
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue<Op<T1, op_scalar_times>, Op<T2, op_scalar_times>, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.m);
  const unwrap<T2> tmp2(in.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const eT k1 = in.A.aux;
  const eT k2 = in.B.aux;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const u32 local_n_elem = A.n_elem;
  
  for(u32 i=0; i<local_n_elem; ++i)
    {
    out_mem[i] = k1*A_mem[i] + k2*B_mem[i];
    }
  
  }



template<typename T1, typename T2, typename T3>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue< Glue<Op<T1, op_scalar_times>, Op<T2, op_scalar_times>, glue_plus>, Op<T3, op_scalar_times>, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.A.m);
  const unwrap<T2> tmp2(in.A.B.m);
  const unwrap<T3> tmp3(in.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  arma_debug_assert_same_size(B, C, "matrix addition");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const eT k1 = in.A.A.aux;
  const eT k2 = in.A.B.aux;
  const eT k3 = in.B.aux;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  const eT* C_mem   = C.mem;
  
  const u32 local_n_elem = A.n_elem;
  
  for(u32 i=0; i<local_n_elem; ++i)
    {
    out_mem[i] = k1*A_mem[i] + k2*B_mem[i] + k3*C_mem[i];
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue<Op<T1, op_scalar_div_pre>, Op<T2, op_scalar_div_pre>, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.m);
  const unwrap<T2> tmp2(in.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const eT k1 = in.A.aux;
  const eT k2 = in.B.aux;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const u32 local_n_elem = A.n_elem;
  
  for(u32 i=0; i<local_n_elem; ++i)
    {
    out_mem[i] = k1/A_mem[i] + k2/B_mem[i];
    }
  
  }



template<typename T1, typename T2, typename T3>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue< Glue<Op<T1, op_scalar_div_pre>, Op<T2, op_scalar_div_pre>, glue_plus>, Op<T3, op_scalar_div_pre>, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.A.m);
  const unwrap<T2> tmp2(in.A.B.m);
  const unwrap<T3> tmp3(in.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  arma_debug_assert_same_size(B, C, "matrix addition");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const eT k1 = in.A.A.aux;
  const eT k2 = in.A.B.aux;
  const eT k3 = in.B.aux;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  const eT* C_mem   = C.mem;
  
  const u32 local_n_elem = A.n_elem;
  
  for(u32 i=0; i<local_n_elem; ++i)
    {
    out_mem[i] = k1/A_mem[i] + k2/B_mem[i] + k3/C_mem[i];
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue<Op<T1, op_scalar_div_post>, Op<T2, op_scalar_div_post>, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.m);
  const unwrap<T2> tmp2(in.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const eT k1 = in.A.aux;
  const eT k2 = in.B.aux;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const u32 local_n_elem = A.n_elem;
  
  for(u32 i=0; i<local_n_elem; ++i)
    {
    out_mem[i] = A_mem[i]/k1 + B_mem[i]/k2;
    }
  
  }



template<typename T1, typename T2, typename T3>
inline
void
glue_plus::apply
  (
  Mat<typename T1::elem_type>& out,
  const Glue< Glue<Op<T1, op_scalar_div_post>, Op<T2, op_scalar_div_post>, glue_plus>, Op<T3, op_scalar_div_post>, glue_plus>& in
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(in.A.A.m);
  const unwrap<T2> tmp2(in.A.B.m);
  const unwrap<T3> tmp3(in.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  const Mat<eT>& C = tmp3.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  arma_debug_assert_same_size(B, C, "matrix addition");
  
  out.set_size(A.n_rows, A.n_cols);
  
  const eT k1 = in.A.A.aux;
  const eT k2 = in.A.B.aux;
  const eT k3 = in.B.aux;
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  const eT* C_mem   = C.mem;
  
  const u32 local_n_elem = A.n_elem;
  
  for(u32 i=0; i<local_n_elem; ++i)
    {
    out_mem[i] = A_mem[i]/k1 + B_mem[i]/k2 + C_mem[i]/k3;
    }
  
  }



//
// matrix addition with different element types

template<typename eT1, typename eT2>
inline
void
glue_plus::apply_mixed(Mat<typename promote_type<eT1,eT2>::result>& out, const Mat<eT1>& X, const Mat<eT2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  arma_debug_assert_same_size(X,Y, "matrix addition");  
  
  out.set_size(X.n_rows, X.n_cols);
  
        out_eT* out_mem = out.memptr();
  const eT1*    X_mem   = X.mem;
  const eT2*    Y_mem   = Y.mem;
  
  const u32 n_elem = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(X_mem[i]) + upgrade_val<eT1,eT2>::apply(Y_mem[i]);
    }
  }



//
// glue_plus_diag


template<typename T1, typename T2>
inline
void
glue_plus_diag::apply(Mat<typename T1::elem_type>& out, const T1& A_orig, const Op<T2,op_diagmat>& B_orig)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(A_orig);
  const unwrap<T2> tmp2(B_orig.m);
  
  const Mat<eT>& A = tmp1.M;  
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( !B.is_square(), "glue_plus_diag::apply(): matrices must be square" );
  arma_debug_assert_same_size(A, B, "matrix addition");

  
  // no aliasing problem
  out.set_size(A.n_rows, A.n_cols);
  
  for(u32 col=0; col<A.n_cols; ++col)
    {
    for(u32 row=0; row<A.n_rows; ++row)
      {
      if(col != row)
        {
        out.at(row,col) = A.at(row,col);
        }
      else
        {
        out.at(row,col) = A.at(row,col) + B.at(row,col);
        }
      }
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_plus_diag::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& A_orig, const Op<T2,op_diagmat>& B_orig)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  const unwrap<T1> tmp1(A_orig.m);
  const unwrap<T2> tmp2(B_orig.m);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( !A.is_square(), "glue_plus_diag::apply(): matrices must be square" );
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  
  if( (&out != &A) && (&out != &B) )
    {
    out.zeros(A.n_rows, A.n_cols);
    
    for(u32 i=0; i<A.n_rows; ++i)
      {
      out.at(i,i) = A.at(i,i) + B.at(i,i);
      }
    }
  else
    {
    out.set_size(A.n_rows, A.n_cols);
  
    for(u32 col=0; col<A.n_cols; ++col)
      {
      for(u32 row=0; row<A.n_rows; ++row)
        {
        if(col != row)
          {
          out.at(row,col) = 0.0;
          }
        else
          {
          out.at(row,col) = A.at(row,col) + B.at(row,col);
          }
        }
      }
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_plus_diag::apply(Mat<typename T1::elem_type>& out, const Glue<T1, Op<T2,op_diagmat>, glue_plus_diag>& X)
  {
  glue_plus_diag::apply(out, X.A, X.B);
  }



template<typename T1, typename T2>
inline
void
glue_plus_diag::apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, T2, glue_plus_diag>& X)
  {
  glue_plus_diag::apply(out, X.B, X.A);  // NOTE: arguments are swapped
  }



template<typename T1, typename T2>
inline
void
glue_plus_diag::apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_plus_diag>& X)
  {
  glue_plus_diag::apply(out, X.A, X.B);
  }



//! @}
