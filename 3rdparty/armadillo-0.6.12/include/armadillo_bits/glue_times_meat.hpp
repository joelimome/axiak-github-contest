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


//! \addtogroup glue_times
//! @{


template<typename eT>
arma_inline
u32 glue_times::mul_storage_cost(const Mat<eT>& X, const Mat<eT>& Y)
  {
  return X.n_rows * Y.n_cols;
  }



//! multiply matrices A and B, storing the result in 'out'
//! assumes that A and B are not aliases of 'out'
template<typename eT>
inline
void
glue_times::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_mul_size(A, B, "matrix multiply");
  
  out.set_size(A.n_rows,B.n_cols);
  gemm<>::apply(out,A,B);
  }



template<typename eT>
inline
void
glue_times::apply(Mat<eT>& out, const Mat<eT>& A_in, const Mat<eT>& B_in)
  {
  arma_extra_debug_sigprint();
  
  if( (&out != &A_in) && (&out != &B_in) )
    {
    glue_times::apply_noalias(out,A_in,B_in);
    }
  else
    {
    
    if( (&out == &A_in) && (&out != &B_in) )
      {
      Mat<eT> A_copy(A_in);
      glue_times::apply_noalias(out,A_copy,B_in);
      }
    else
    if( (&out != &A_in) && (&out == &B_in) )
      {
      Mat<eT> B_copy(B_in);
      glue_times::apply_noalias(out,A_in,B_copy);
      }
    else
    if( (&out == &A_in) && (&out == &B_in) )
      {
      Mat<eT> tmp(A_in);
      glue_times::apply_noalias(out,tmp,tmp);
      }

    }
    
  }


template<typename eT>
inline
void
glue_times::apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C)
  {
  arma_extra_debug_sigprint();

  arma_debug_assert_mul_size(A, B, "matrix multiply");
  arma_debug_assert_mul_size(B, C, "matrix multiply");
  
  if( mul_storage_cost(A,B) <= mul_storage_cost(B,C) )
    {
    Mat<eT> tmp;
    glue_times::apply_noalias(tmp, A, B);
    
    if(&out != &C)
      {
      glue_times::apply_noalias(out, tmp, C);
      }
    else
      {
      Mat<eT> C_copy = C;
      glue_times::apply_noalias(out, tmp, C_copy);
      }
      
    }
  else
    {
    Mat<eT> tmp;
    glue_times::apply_noalias(tmp, B, C);
    
    if(&out != &A)
      {
      glue_times::apply_noalias(out, A, tmp);
      }
    else
      {
      Mat<eT> A_copy = A;
      glue_times::apply_noalias(out, A_copy, tmp);
      }
    }
  
  }



template<typename T1, typename T2>
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const s32 N_mat = 1 + depth_lhs< glue_times, Glue<T1,T2,glue_times> >::num;

  arma_extra_debug_print(arma_boost::format("N_mat = %d") % N_mat);

  if(N_mat == 2)
    {
    const unwrap<T1> tmp1(X.A);
    const unwrap<T2> tmp2(X.B);
    
    glue_times::apply(out, tmp1.M, tmp2.M);
    }
  else
    {
    // we have at least three matrices

    const Mat<eT>* ptrs[N_mat];
    bool            del[N_mat];
  
    // takes care of any aliasing problems
    mat_ptrs_outcheck<glue_times, Glue<T1,T2,glue_times> >::get_ptrs(ptrs, del, X, &out);
  
    for(s32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format("ptrs[%d] = %x") % i % ptrs[i] );
    for(s32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format(" del[%d] = %d") % i %  del[i] );
  
  
    arma_extra_debug_print( arma_boost::format("required size of 'out':  %d, %d") % ptrs[0]->n_rows % ptrs[N_mat-1]->n_cols );
  
    int order[N_mat];  for(s32 i=0; i<N_mat; ++i)  order[i] = -1;
  
    int first_id = 0;
    int last_id  = N_mat-1;
    int starting_id = -1;
  
    int mat_count = N_mat;
  
    int largest_size = 0;
  
    while(mat_count != 0)
      {
  
      for(s32 i=first_id; i != N_mat; ++i)
        {
        if(order[i] == -1)  { first_id = i; break; }
        }
  
      for(s32 i=last_id; i != -1; --i)
        {
        if(order[i] == -1)  { last_id = i; break; }
        }
  
      arma_extra_debug_print();
      arma_extra_debug_print(arma_boost::format("mat_count = %d") % mat_count );
      arma_extra_debug_print(arma_boost::format("first_id  = %d") % first_id  );
      arma_extra_debug_print(arma_boost::format("last_id   = %d") % last_id   );
  
      if(first_id == last_id)  { order[first_id] = 0; starting_id = first_id; break; }
  
      s32 storage_cost_wo_last  = mul_storage_cost( *ptrs[ first_id   ], *ptrs[ last_id-1 ] );
      s32 storage_cost_wo_first = mul_storage_cost( *ptrs[ first_id+1 ], *ptrs[ last_id   ] );
  
      if(storage_cost_wo_last < storage_cost_wo_first)
        {
        order[last_id]  = mat_count-1;
        if(storage_cost_wo_last > largest_size)  largest_size = storage_cost_wo_last;
        }
      else
        {
        order[first_id] = mat_count-1;
        if(storage_cost_wo_first > largest_size)  largest_size = storage_cost_wo_first;
        }
  
      arma_extra_debug_print(arma_boost::format("storage_cost_wo_last  = %d") % storage_cost_wo_last  );
      arma_extra_debug_print(arma_boost::format("storage_cost_wo_first = %d") % storage_cost_wo_first );
  
      arma_extra_debug_print("order = ");
      for(s32 i=0; i != N_mat; ++i)  arma_extra_debug_print(order[i]);
  
      --mat_count;
      }
  
    arma_extra_debug_print("final order = ");
    for(s32 i=0; i != N_mat; ++i)  arma_extra_debug_print(order[i]);
  
    arma_extra_debug_print(arma_boost::format("*** largest_size = %d") % largest_size);
    arma_extra_debug_print(arma_boost::format("starting_id  = %d") % starting_id);
  
  
    // multiply based on order
    // if there are only three matrices, we need only one temporary store:
    //   out = a*b*c translates to:  tmp1 = a*b,  out = tmp1*c
    //
    // if there are four matrices, we need two temporary stores
    //   out = a*b*c*d translates to:  tmp1 = a*b, tmp2 = tmp1*c, out = tmp2*d
    //
    // if there are five matrices, we need two temporary stores
    //   out = a*b*c*d*e translates to:  tmp1 = a*b, tmp2 = tmp1*c, tmp1 = tmp2*d, out = tmp1*e
    //
    // if there are six matrices, we need two temporary stores
    //   out = a*b*c*d*e*f translates to:  tmp1 = a*b, tmp2 = tmp1*c, tmp1 = tmp2*d, tmp2 = tmp1*e, out = tmp2*f
    //
  
    
    const u32 N_mul = N_mat - 1;
    int mul_count = N_mul;
    int current_id = starting_id;
  
    const Mat<eT>* src_mat_1_ptr = ptrs[current_id];
    const Mat<eT>* src_mat_2_ptr = 0;
  
    // TODO:
    // allocate two storage areas (of size 'largest_size'), not two matrices
  
    
    Mat<eT> tmp_mat_1;
    Mat<eT> tmp_mat_2;
    
    Mat<eT>* tmp_mat_1_ptr = &tmp_mat_1;
    Mat<eT>* tmp_mat_2_ptr = (N_mul <= 2) ? 0 : &tmp_mat_2;
  
    Mat<eT>* dest_mat_ptr  = tmp_mat_2_ptr;
  
    arma_extra_debug_print(arma_boost::format("tmp_mat_1_ptr = %x") % tmp_mat_1_ptr );
    arma_extra_debug_print(arma_boost::format("tmp_mat_2_ptr = %x") % tmp_mat_2_ptr );
    arma_extra_debug_print(arma_boost::format("&out          = %x") % &out );
  
    while(mul_count != 0)
      {
      arma_extra_debug_print("");
      arma_extra_debug_print("");
      arma_extra_debug_print(arma_boost::format("mul_count = %d") % mul_count);
  
      arma_extra_debug_print("order = ");
      for(s32 i=0; i != N_mat; ++i)  arma_extra_debug_print(order[i]);
      arma_extra_debug_print("");
  
      // only one multiplication left, hence destination matrix is the out matrix
      if(mul_count == 1)
        {
        arma_extra_debug_print("dest_mat = &out");
        dest_mat_ptr = &out;
        }
      else
        {
        if(dest_mat_ptr == tmp_mat_2_ptr)
          {
          arma_extra_debug_print("dest_mat_ptr = tmp_mat_2_ptr");
          dest_mat_ptr = tmp_mat_1_ptr;
          }
        else
          {
          arma_extra_debug_print("dest_mat_ptr = tmp_mat_1_ptr");
          dest_mat_ptr = tmp_mat_2_ptr;
          }
        }
  
      arma_extra_debug_print(arma_boost::format("dest_mat_ptr = %x") % dest_mat_ptr );
  
      // search on either side of current_pos for a useable value.  unuseable values are equal to -1
      s32 left_val = N_mat;
      s32 left_id = -1;
  
      s32 right_val = N_mat;
      s32 right_id = -1;
  
      // go left from current_pos
      for(s32 i=current_id-1; i >= 0; --i)
        if( order[i] > order[current_id] ) { left_val = order[i]; left_id = i; break; }
  
      // go right from current_pos
      for(s32 i=current_id+1; i < N_mat; ++i)
        if( order[current_id] < order[i] ) { right_val = order[i]; right_id = i; break; }
  
      arma_extra_debug_print("");
      arma_extra_debug_print(arma_boost::format("left_id  = %d") % left_id  );
      arma_extra_debug_print(arma_boost::format("left_val = %f") % left_val );
  
      arma_extra_debug_print("");
      arma_extra_debug_print(arma_boost::format("right_id  = %d") % right_id  );
      arma_extra_debug_print(arma_boost::format("right_val = %f") % right_val );
  
  
      if(left_val < right_val)
        {
        // a pre-multiply
        src_mat_2_ptr = ptrs[left_id];
  
        arma_extra_debug_print("");
        arma_extra_debug_print(arma_boost::format("case pre-multiply with matrix %d") % left_id);
        arma_extra_debug_print(arma_boost::format("required destination size: %d, %d  (%d)") %  src_mat_2_ptr->n_rows % src_mat_1_ptr->n_cols % (src_mat_2_ptr->n_rows * src_mat_1_ptr->n_cols) );
  
        glue_times::apply_noalias(*dest_mat_ptr, *src_mat_2_ptr, *src_mat_1_ptr);
  
        order[current_id] = -1;
        current_id = left_id;
        }
      else
        {
        // a post-multiply
        src_mat_2_ptr = ptrs[right_id];
  
        arma_extra_debug_print("");
        arma_extra_debug_print(arma_boost::format("case post-multiply with matrix %d") % right_id);
        arma_extra_debug_print(arma_boost::format("required destination size: %d, %d  (%d)") % src_mat_1_ptr->n_rows % src_mat_2_ptr->n_cols % (src_mat_1_ptr->n_rows * src_mat_2_ptr->n_cols) );
  
        glue_times::apply_noalias(*dest_mat_ptr, *src_mat_1_ptr, *src_mat_2_ptr);
  
        order[current_id] = -1;
        current_id = right_id;
        }
  
      // update pointer to source matrix: must point to last multiplication result
      src_mat_1_ptr = dest_mat_ptr;
  
      --mul_count;
      }
  
  
    for(s32 i=0; i<N_mat; ++i)
      {
      if(del[i] == true)
        {
        arma_extra_debug_print(arma_boost::format("delete mat_ptr[%d]") % i );
        delete ptrs[i];
        }
      }
    }
  }



template<typename eT>
inline
void
glue_times::apply(Mat<eT>& out, const Glue<Mat<eT>,Mat<eT>,glue_times>& X)
  {
  glue_times::apply(out, X.A, X.B);
  }



template<typename eT>
inline
void
glue_times::apply(Mat<eT>& out, const Glue< Glue<Mat<eT>,Mat<eT>, glue_times>, Mat<eT>, glue_times>& X)
  {
  glue_times::apply(out, X.A.A, X.A.B, X.B);
  }



template<typename eT>
inline
void
glue_times::apply_inplace(Mat<eT>& out, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_mul_size(out, B, "matrix multiply");
  
  if(out.n_cols == B.n_cols)
    {
    podarray<eT> tmp(out.n_cols);
    eT* tmp_rowdata = tmp.memptr();
    
    for(u32 out_row=0; out_row < out.n_rows; ++out_row)
      {
      for(u32 out_col=0; out_col < out.n_cols; ++out_col)
        {
        tmp_rowdata[out_col] = out.at(out_row,out_col);
        }
      
      for(u32 B_col=0; B_col < B.n_cols; ++B_col)
        {
        const eT* B_coldata = B.colptr(B_col);
        
        eT val = eT(0);
        for(u32 i=0; i < B.n_rows; ++i)
          {
          val += tmp_rowdata[i] * B_coldata[i];
          }
        
        out.at(out_row,B_col) = val;
        }
      }
    
    }
  else
    {
    Mat<eT> tmp = out;
    glue_times::apply(out, tmp, B);
    }
  
  }



template<typename T1, typename op_type>
inline
void
glue_times::apply_inplace(Mat<typename T1::elem_type>& out, const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT> tmp(X);
  glue_times::apply(out, out, tmp);
  }



template<typename T1, typename T2, typename glue_type>
inline
void
glue_times::apply_inplace(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  out = out * X;
  }



//! out = T1 * trans(T2)
template<typename T1, typename T2>
inline
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue<T1, Op<T2,op_trans>, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  
  typedef typename T1::elem_type eT;
  
  // checks for aliases are done later
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B.m);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_cols, B.n_rows, "matrix multiply");
    
  if( (A.n_rows*B.n_rows) > 0)
    {
    if(&A != &B)   // A*B'
      {
      unwrap_check< Mat<eT> > A_safe_tmp(A, out);
      unwrap_check< Mat<eT> > B_safe_tmp(B, out);
      
      const Mat<eT>& A_safe = A_safe_tmp.M;
      const Mat<eT>& B_safe = B_safe_tmp.M;
      
      out.set_size(A_safe.n_rows, B_safe.n_rows);
      
      gemm<false,true>::apply(out, A, B);
      }
    else   // A*A'
      {
      arma_extra_debug_print("glue_times::apply(): detected A*A'");
      
      Mat<eT> tmp;
      op_trans::apply(tmp,A);
      
      // no aliasing problem
      out.set_size(A.n_rows, A.n_rows);
      
      for(u32 row=0; row != A.n_rows; ++row)
        {
        for(u32 col=0; col <= row; ++col)
          {
          const eT* coldata1 = tmp.colptr(row);
          const eT* coldata2 = tmp.colptr(col);
        
          eT val = eT(0);
          for(u32 i=0; i < tmp.n_rows; ++i)
            {
            val += coldata1[i] * coldata2[i];
            }
          
          out.at(row,col) = val;
          out.at(col,row) = val;
          }
        }
        
      }
    
    }
  
  }



//! out = trans(T1) * T2
template<typename T1, typename T2>
inline
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue< Op<T1,op_trans>, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp1(X.A.m, out);
  const unwrap_check<T2> tmp2(X.B,   out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A.n_cols, A.n_rows, B.n_rows, B.n_cols, "matrix multiply");
  
  if( (A.n_cols*B.n_cols) > 0 )
    {
    out.set_size(A.n_cols, B.n_cols);
    
    gemm<true,false>::apply(out, A, B);
    }
    
  }



//! out = trans(T1) * trans(T2)
template<typename T1, typename T2>
inline
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue< Op<T1,op_trans>, Op<T2,op_trans>, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp1(X.A.m, out);
  const unwrap_check<T2> tmp2(X.B.m, out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A.n_cols, A.n_rows, B.n_cols, B.n_rows, "matrix multiply");
  
  if( (A.n_cols*B.n_rows) > 0 )
    {
    out.set_size(A.n_cols, B.n_rows);
    
    gemm<true,true>::apply(out, A, B);
    
    }
    
  }




//! out = -T1 * T2
template<typename T1, typename T2>
inline
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue< Op<T1, op_neg>, T2, glue_times>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> tmp1(X.A.m, out);
  const unwrap_check<T2> tmp2(X.B,   out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;

  glue_times::apply(out, A, B);
  
  const u32 n_elem = out.n_elem;
  for(u32 i=0; i<n_elem; ++i)
    {
    out[i] = -out[i];
    }
  }




template<typename eT>
inline
eT
glue_times::direct_rowvec_mat_colvec
  (
  const eT*            A_mem,
  const Mat<eT>& B,
  const eT*            C_mem
  )
  {
  arma_extra_debug_sigprint();
  
  const u32 cost_AB = B.n_cols;
  const u32 cost_BC = B.n_rows;
  
  if(cost_AB <= cost_BC)
    {
    podarray<eT> tmp(B.n_cols);
    
    for(u32 col=0; col<B.n_cols; ++col)
      {
      const eT* B_coldata = B.colptr(col);
      
      eT val = eT(0);
      for(u32 i=0; i<B.n_rows; ++i)
        {
        val += A_mem[i] * B_coldata[i];
        }
        
      tmp[col] = val;
      }
    
    return op_dot::direct_dot(B.n_cols, tmp.mem, C_mem);
    }
  else
    {
    podarray<eT> tmp(B.n_rows);
    
    for(u32 row=0; row<B.n_rows; ++row)
      {
      eT val = eT(0);
      for(u32 col=0; col<B.n_cols; ++col)
        {
        val += B.at(row,col) * C_mem[col];
        }
      
      tmp[row] = val;
      }
    
    return op_dot::direct_dot(B.n_rows, A_mem, tmp.mem);
    }
  
  
  }



template<typename eT>
inline
eT
glue_times::direct_rowvec_diagmat_colvec
  (
  const eT*            A_mem,
  const Mat<eT>& B,
  const eT*            C_mem
  )
  {
  arma_extra_debug_sigprint();
  
  eT val = eT(0);

  for(u32 i=0; i<B.n_rows; ++i)
    {
    val += A_mem[i] * B.at(i,i) * C_mem[i];
    }

  return val;
  }



template<typename eT>
inline
eT
glue_times::direct_rowvec_invdiagmat_colvec
  (
  const eT*            A_mem,
  const Mat<eT>& B,
  const eT*            C_mem
  )
  {
  arma_extra_debug_sigprint();
  
  eT val = eT(0);

  for(u32 i=0; i<B.n_rows; ++i)
    {
    val += (A_mem[i] * C_mem[i]) / B.at(i,i);
    }

  return val;
  }



template<typename eT>
inline
eT
glue_times::direct_rowvec_invdiagvec_colvec
  (
  const eT*            A_mem,
  const Mat<eT>& B,
  const eT*            C_mem
  )
  {
  arma_extra_debug_sigprint();
  
  const eT* B_mem = B.mem;
  
  eT val = eT(0);

  for(u32 i=0; i<B.n_elem; ++i)
    {
    val += (A_mem[i] * C_mem[i]) / B_mem[i];
    }

  return val;
  }



//
// matrix multiplication with different element types

template<typename eT1, typename eT2>
inline
void
glue_times::apply_mixed(Mat<typename promote_type<eT1,eT2>::result>& out, const Mat<eT1>& X, const Mat<eT2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  arma_debug_assert_mul_size(X,Y, "matrix multiply");
  
  out.set_size(X.n_rows,Y.n_cols);
  gemm_mixed<>::apply(out, X, Y);
  }



//
// glue_times_diag


template<typename T1, typename T2>
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const T1& A_orig, const Op<T2,op_diagmat>& B_orig)
  {
  arma_extra_debug_sigprint();
    
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  const unwrap_check<T1> tmp1(A_orig,   out);
  const unwrap_check<T2> tmp2(B_orig.m, out);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (B.is_square() == false), "glue_times_diag::apply(): incompatible matrix dimensions" );
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiply");
  
  out.set_size(A.n_rows, B.n_cols);
  
  for(u32 col=0; col<A.n_cols; ++col)
    {
    const eT val = B.at(col,col);
    
    const eT*   A_coldata =   A.colptr(col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row=0; row<B.n_rows; ++row)
      {
      out_coldata[row] = A_coldata[row] * val;
      }
    
    }
  
  }



template<typename T1, typename T2>
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& A_orig, const T2& B_orig)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  const unwrap_check<T1> tmp1(A_orig.m, out);
  const unwrap_check<T2> tmp2(B_orig,   out);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( (A.is_square() == false), "glue_times_diag::apply(): incompatible matrix dimensions" );
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiply");
  
  out.set_size(A.n_rows, B.n_cols);
  
  
  for(u32 col=0; col<A.n_cols; ++col)
    {
    const eT*   B_coldata =   B.colptr(col);
          eT* out_coldata = out.colptr(col);
    
    for(u32 row=0; row<B.n_rows; ++row)
      {
      out_coldata[row] = A.at(row,row) * B_coldata[row];
      }
    
    }

  }



template<typename T1, typename T2>
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& A_orig, const Op<T2,op_diagmat>& B_orig)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
  
  unwrap_check<T1> tmp1(A_orig.m, out);
  unwrap_check<T2> tmp2(B_orig.m, out);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check( !A.is_square() || !B.is_square(), "glue_times_diag::apply(): incompatible matrix dimensions" );
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiply");
  
  out.zeros(A.n_rows, B.n_cols);
  
  for(u32 i=0; i<A.n_rows; ++i)
    {
    out.at(i,i) = A.at(i,i) * B.at(i,i);
    }
  }



template<typename T1, typename T2>
inline
void 
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const Glue<T1, Op<T2,op_diagmat>, glue_times_diag>& X)
  {
  glue_times_diag::apply(out, X.A, X.B);
  }



template<typename T1, typename T2>
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, T2, glue_times_diag>& X)
  {
  glue_times_diag::apply(out, X.A, X.B);
  }



template<typename T1, typename T2>
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_times_diag>& X)
  {
  glue_times_diag::apply(out, X.A, X.B);
  }



//
// glue_times_vec



template<typename eT>
inline
void
glue_times_vec::mul_col_row(Mat<eT>& out, const eT* A, const eT* B)
  {
  const u32 n_rows = out.n_rows;
  const u32 n_cols = out.n_cols;
  
  for(u32 col=0; col < n_cols; ++col)
    {
    const eT val = B[col];
    
    eT* out_coldata = out.colptr(col);
    
    for(u32 row=0; row < n_rows; ++row)
      {
      out_coldata[row] = A[row] * val;
      }
    }
  
  }



template<typename T1>
inline
void
glue_times_vec::apply(Mat<typename T1::elem_type>& out, const Glue<T1, Col<typename T1::elem_type>,glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  unwrap_check< T1 >               tmp1(X.A, out);
  unwrap_check< Col<eT> > tmp2(X.B, out);
  
  const Mat<eT>& A = tmp1.M;
  const Col<eT>& B = tmp2.M;

  arma_debug_assert_mul_size(A, B, "vector multiply");
  
  out.set_size(A.n_rows, 1);
  
  //gemm<>::apply(out,A,B);  // NOTE: B is interpreted as a Mat
  gemv<>::apply(out.memptr(), A, B.mem);
  }



template<typename T1>
inline
void
glue_times_vec::apply(Mat<typename T1::elem_type>& out, const Glue<T1, Row<typename T1::elem_type>,glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  // T1 * rowvec makes sense only if T1 ends up being a matrix with one column (i.e. a column vector)
  
  typedef typename T1::elem_type eT;
  
  unwrap_check< T1 >               tmp1(X.A, out);
  unwrap_check< Row<eT> > tmp2(X.B, out);
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;   // NOTE: interpretation of a Row as a Mat
  
  arma_debug_assert_mul_size(A, B, "vector multiply");
  
  out.set_size(A.n_rows, B.n_cols);
  
  glue_times_vec::mul_col_row(out, A.mem, B.mem);
  }



template<typename T1>
inline
void
glue_times_vec::apply(Mat<typename T1::elem_type>& out, const Glue<Col<typename T1::elem_type>,T1,glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  // colvec * T1 makes sense only if T1 ends up being a matrix with one row (i.e. a row vector)
  
  
  typedef typename T1::elem_type eT;
  
  unwrap_check< Col<eT> > tmp1(X.A, out);
  unwrap_check< T1 >               tmp2(X.B, out);
  
  const Mat<eT>& A = tmp1.M;   // NOTE: interpretation of a Col as a Mat
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "vector multiply");
  
  out.set_size(A.n_rows, B.n_cols);
  
  glue_times_vec::mul_col_row(out, A.mem, B.mem);
  }



template<typename T1>
inline
void
glue_times_vec::apply(Mat<typename T1::elem_type>& out, const Glue<Row<typename T1::elem_type>,T1,glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  unwrap_check< Row<eT> > tmp1(X.A, out);
  unwrap_check< T1 >               tmp2(X.B, out);
  
  const Row<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "vector multiply");
  
  out.set_size(A.n_rows, B.n_cols);
  
//         eT* out_mem = out.memptr();
//   const eT* A_mem   = A.mem;
//   
//   const u32 A_n_cols = A.n_cols;
//   const u32 B_n_cols = B.n_cols;
//   
//   for(u32 col=0; col<B_n_cols; ++col)
//     {
//     const eT* B_coldata = B.colptr(col);
//     
//     eT val = eT(0);
//     for(u32 i=0; i<A_n_cols; ++i)
//       {
//       val += A_mem[i] * B_coldata[i];
//       }
//       
//     out_mem[col] = val;
//     }
  
  gemv<true>::apply(out.memptr(), B, A.mem);
  }



template<typename eT>
inline
void
glue_times_vec::apply(Mat<eT>& out, const Glue<Col<eT>,Row<eT>,glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  unwrap_check< Col<eT> > tmp1(X.A, out);
  unwrap_check< Row<eT> > tmp2(X.B, out);
  
  const Col<eT>& A = tmp1.M;
  const Row<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A, B, "vector multiply");
  
  out.set_size(A.n_rows, B.n_cols);
  
  glue_times_vec::mul_col_row(out, A.mem, B.mem);
  }



template<typename eT>
inline
void
glue_times_vec::apply(Mat<eT>& out, const Glue< Op<Row<eT>, op_trans>, Row<eT>, glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  unwrap_check< Row<eT> > tmp1(X.A.m, out);
  unwrap_check< Row<eT> > tmp2(X.B,   out);
  
  const Row<eT>& A = tmp1.M;
  const Row<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A.n_cols, A.n_rows, B.n_rows, B.n_cols, "vector multiply");
  
  out.set_size(A.n_cols, B.n_cols);
  
  glue_times_vec::mul_col_row(out, A.mem, B.mem);
  }



template<typename eT>
inline
void
glue_times_vec::apply(Mat<eT>& out, const Glue< Col<eT>, Op<Col<eT>, op_trans>, glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  unwrap_check< Col<eT> > tmp1(X.A,   out);
  unwrap_check< Col<eT> > tmp2(X.B.m, out);
  
  const Col<eT>& A = tmp1.M;
  const Col<eT>& B = tmp2.M;
  
  arma_debug_assert_mul_size(A.n_rows, A.n_cols, B.n_cols, B.n_rows, "vector multiply");
  
  out.set_size(A.n_rows, B.n_rows);
  
  glue_times_vec::mul_col_row(out, A.mem, B.mem);
  }



template<typename T1>
inline
void
glue_times_vec::apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1, op_trans>, Col<typename T1::elem_type>,glue_times_vec>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  unwrap_check< T1 >               tmp1(X.A.m, out);
  unwrap_check< Col<eT> > tmp2(X.B,   out);
  
  const Mat<eT>& A = tmp1.M;
  const Col<eT>& B = tmp2.M;

  arma_debug_assert_mul_size(A.n_cols, A.n_rows, B.n_rows, B.n_cols, "vector multiply");
  
  out.set_size(A.n_cols, B.n_cols);
  
//         eT* out_mem = out.memptr();
//   const eT* B_mem   = B.mem;
//   
//   const u32 A_n_cols = A.n_cols;
//   const u32 B_n_rows = B.n_rows;
//   
//   for(u32 col=0; col < A_n_cols; ++col)
//     {
//     const eT* A_col = A.colptr(col);
//     
//     eT val = eT(0);
//     for(u32 row=0; row<B_n_rows; ++row)
//       {
//       val += A_col[row] * B_mem[row];
//       }
//     
//     out_mem[col] = val;
//     }
  
  gemv<true>::apply(out.memptr(), A, B.mem);
  }



//! @}
