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


//! \addtogroup glue_div
//! @{


template<typename eT>
inline
void
glue_div::apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "element-wise matrix division");
    
  // no aliasing problem
  out.set_size(A.n_rows, A.n_cols);
    
  const u32 n_elem = A.n_elem;
    
  for(u32 i=0; i<n_elem; ++i)
    {
    out[i] = A.mem[i] / B.mem[i];
    }
  
  }



template<typename eT>
inline
void
glue_div::apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(A, B, "element-wise matrix division");
  arma_debug_assert_same_size(A, C, "element-wise matrix division");
  
  // no aliasing problem
  out.set_size(A.n_rows, A.n_cols);
    
  const u32 n_elem = A.n_elem;
  for(u32 i=0; i != n_elem; ++i)
    {
    out[i] = A.mem[i] / B.mem[i] / C.mem[i];
    }
  
  }



template<typename eT>
inline
void
glue_div::apply(Mat<eT>& out, const Glue<Mat<eT>,Mat<eT>,glue_div>& X)
  {
  glue_div::apply(out, X.A, X.B);
  }



template<typename eT>
inline
void
glue_div::apply(Mat<eT>& out, const Glue< Glue<Mat<eT>,Mat<eT>,glue_div>, Mat<eT>,glue_div>& X)
  {
  glue_div::apply(out, X.A.A, X.A.B, X.B);
  }



template<typename T1, typename T2>
void
glue_div::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_div>& X)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const u32 N_mat = 1 + depth_lhs< glue_div, Glue<T1,T2,glue_div> >::num;
  arma_extra_debug_print( arma_boost::format("N_mat = %d") % N_mat );

  if(N_mat == 2)
    {
    const unwrap<T1> tmp1(X.A);
    const unwrap<T2> tmp2(X.B);
    
    glue_div::apply(out, tmp1.M, tmp2.M);
    }
  else
    {
    const Mat<eT>* ptrs[N_mat];
    bool            del[N_mat];

    mat_ptrs<glue_div, Glue<T1,T2,glue_div> >::get_ptrs(ptrs, del, X);
    //mat_ptrs_outcheck<glue_div, Glue<T1,T2,glue_div> >::get_ptrs(ptrs, del, X, &out);

    for(u32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format("ptrs[%d] = %x") % i % ptrs[i] );
    for(u32 i=0; i<N_mat; ++i)  arma_extra_debug_print( arma_boost::format(" del[%d] = %d") % i %  del[i] );

    const Mat<eT>& tmp_mat = *(ptrs[0]);

    for(u32 i=1; i<N_mat; ++i)
      {
      arma_debug_assert_same_size(tmp_mat, *(ptrs[i]), "element-wise matrix division");
      }
  
    
    const u32 n_rows = ptrs[0]->n_rows;
    const u32 n_cols = ptrs[0]->n_cols;

    // no aliasing problem
    out.set_size(n_rows,n_cols);
    
    const u32 n_elem = ptrs[0]->n_elem;
    
    for(u32 j=0; j<n_elem; ++j)
      {
      eT acc = ptrs[0]->mem[j];
      
      for(u32 i=1; i<N_mat; ++i)
        {
        acc /= ptrs[i]->mem[j];
        }
      
      out[j] = acc;
      }
    
    
    for(u32 i=0; i<N_mat; ++i)
      {
      if(del[i] == true)
        {
        arma_extra_debug_print( arma_boost::format("delete ptrs[%d]") % i );
        delete ptrs[i];
        }
      }

    }
  }



template<typename eT>
inline
void
glue_div::apply_inplace(Mat<eT>& out, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out, B, "element-wise matrix division");
  
  const u32 n_elem = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out[i] /= B.mem[i];
    }
  
  }



template<typename T1, typename op_type>
inline
void
glue_div::apply_inplace(Mat<typename T1::elem_type>& out, const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT> tmp(X);
  glue_div::apply(out, out, tmp);
  }
  


template<typename T1, typename T2, typename glue_type>
inline
void
glue_div::apply_inplace(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT> tmp(X);
  glue_div::apply(out, X, out);
  }



//
// element-wise division with different element types

template<typename eT1, typename eT2>
inline
void
glue_div::apply_mixed(Mat<typename promote_type<eT1,eT2>::result>& out, const Mat<eT1>& X, const Mat<eT2>& Y)
  {
  arma_extra_debug_sigprint();
  
  typedef typename promote_type<eT1,eT2>::result out_eT;
  
  arma_debug_assert_same_size(X,Y, "element-wise matrix division");
  
  out.set_size(X.n_rows, X.n_cols);
  
        out_eT* out_mem = out.memptr();
  const eT1*    X_mem   = X.mem;
  const eT2*    Y_mem   = Y.mem;
  
  const u32 n_elem = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = upgrade_val<eT1,eT2>::apply(X_mem[i]) / upgrade_val<eT1,eT2>::apply(Y_mem[i]);
    }
  }



//! @}
