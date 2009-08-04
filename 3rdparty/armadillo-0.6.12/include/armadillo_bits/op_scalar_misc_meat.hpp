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



//! \addtogroup op_scalar_misc
//! @{

//! Add a scalar to all elements of a matrix and store the result in a dense matrix
template<typename T1>
inline
void
op_scalar_plus::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_plus>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // no alias problems
  out.set_size(X.n_rows, X.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.mem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out_mem[i] = X_mem[i] + k;
    }

  }



//! For each element of a matrix, subtract it from a scalar and store the result in a dense matrix
template<typename T1>
inline
void
op_scalar_minus_pre::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_minus_pre>& in)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // no alias problems
  out.set_size(X.n_rows, X.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.mem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out_mem[i] = k - X_mem[i];
    }

  }



//! subtract a scalar from each element of a matrix
template<typename T1>
inline
void
op_scalar_minus_post::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_minus_post>& in)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;

  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // no alias problems
  out.set_size(X.n_rows, X.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.mem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out_mem[i] = X_mem[i] - k;
    }

  }



//! Multiply all elements of a matrix by a scalar and store the result in a dense matrix
template<typename T1>
inline
void
op_scalar_times::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // no alias problems
  out.set_size(X.n_rows, X.n_cols);
        
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.mem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i < X.n_elem; ++i)
    {
    out_mem[i] = X_mem[i] * k;
    }
  
  }


  
//! \brief 
//! Evaluate Glue<T1,T2,glue_type>, and then multiply each element of the result by a scalar.
//! Store the final result in a dense matrix.
template<typename T1, typename T2, typename glue_type>
inline
void
op_scalar_times::apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_type>, op_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = in.m;  // implicit conversion to 'Mat<eT>'
  
        eT* out_mem = out.memptr();
  const eT        k = in.aux;
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] *= k;
    }
  }
  
  

//! \brief 
//! Evaluate A + B, and then multiply each element of the result by a scalar.
//! Store the final result in a dense matrix.
template<typename T1, typename T2>
inline
void
op_scalar_times::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_plus>, op_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A_mem[i] + B_mem[i]) * k;
    }
  
  }



//! \brief 
//! Evaluate A - B, and then multiply each element of the result by a scalar.
//! Store the final result in a dense matrix.
template<typename T1, typename T2>
inline
void
op_scalar_times::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_minus>, op_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix subtraction");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A_mem[i] - B_mem[i]) * k;
    }
  
  }



//! \brief 
//! Evaluate A % B (where % is the element-wise multiply operation) and then multiply each element of the result by a scalar.
//! Store the final result in a dense matrix.
template<typename T1, typename T2>
inline
void
op_scalar_times::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_schur>, op_scalar_times>& in)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(in.m.A);
  const unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "schur product");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A_mem[i] * B_mem[i]) * k;
    }
  
  }



//
// 
// 



template<typename T1>
inline
void
op_scalar_div_pre::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_div_pre>& in)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // TODO: analyse effects of aliasing
  out.set_size(X.n_rows, X.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.mem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i < X.n_elem; ++i)
    {
    out_mem[i] = k / X_mem[i];
    }
  
  }


  
template<typename eT>
inline
void
op_scalar_div_pre::apply(Mat<eT>& out, const Op<Mat<eT>,op_scalar_div_pre>& in)
  {
  arma_extra_debug_sigprint();

  const Mat<eT>& X = in.m;
  
  if(&out != &X)
    {
    out.set_size(X.n_rows, X.n_cols);
    
          eT* out_mem = out.memptr();
    const eT* X_mem   = X.mem;
    const eT  k       = in.aux;
    
    for(u32 i=0; i < X.n_elem; ++i)
      {
      out_mem[i] = k / X_mem[i];
      }
    }
  else
    {
          eT* out_mem = out.memptr();
    const eT  k       = in.aux;
    
    for(u32 i=0; i < out.n_elem; ++i)
      {
      out_mem[i] = k / out_mem[i];
      }
    
    }
  
  }



template<typename T1, typename T2, typename glue_type>
inline
void
op_scalar_div_pre::apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_type>, op_scalar_div_pre>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = in.m;  // implicit conversion to 'Mat<eT>'
  
        eT* out_mem = out.memptr();
  const eT  k       = in.aux;
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / out_mem[i];
    }
  }



template<typename T1, typename T2>
inline
void
op_scalar_div_pre::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_plus>, op_scalar_div_pre>& in)
  {
  arma_extra_debug_sigprint();
  
  unwrap<T1> tmp1(in.m.A);
  unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / (A_mem[i] + B_mem[i]);
    }
  
  }



template<typename T1, typename T2>
inline
void
op_scalar_div_pre::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_minus>, op_scalar_div_pre>& in)
  {
  arma_extra_debug_sigprint();
  
  unwrap<T1> tmp1(in.m.A);
  unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix subtraction");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / (A_mem[i] - B_mem[i]);
    }
  
  }



template<typename T1, typename T2>
inline
void
op_scalar_div_pre::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_schur>, op_scalar_div_pre>& in)
  {
  arma_extra_debug_sigprint();
  
  unwrap<T1> tmp1(in.m.A);
  unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "schur product");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = k / (A_mem[i] * B_mem[i]);
    }
  
  }



//
//
//



template<typename T1>
inline
void
op_scalar_div_post::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_div_post>& in)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  const Mat<eT>& X = tmp.M;
  
  // TODO: analyse effects of aliasing
  out.set_size(X.n_rows, X.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* X_mem   = X.mem;
  const eT  k       = in.aux;
  
  for(u32 i=0; i < X.n_elem; ++i)
    {
    out_mem[i] = X_mem[i] / k;
    }
  
  }



template<typename eT>
inline
void
op_scalar_div_post::apply(Mat<eT>& out, const Op<Mat<eT>,op_scalar_div_post>& in)
  {
  arma_extra_debug_sigprint();

  const Mat<eT>& X = in.m;
  
  if(&out != &X)
    {
    out.set_size(X.n_rows, X.n_cols);
    
          eT* out_mem = out.memptr();
    const eT* X_mem   = X.mem;
    const eT  k       = in.aux;
    
    for(u32 i=0; i < X.n_elem; ++i)
      {
      out_mem[i] = X_mem[i] / k;
      }
    }
  else
    {
          eT* out_mem = out.memptr();
    const eT  k       = in.aux;
    
    for(u32 i=0; i < out.n_elem; ++i)
      {
      out_mem[i] /= k;
      }
    
    }
  
  }



template<typename T1, typename T2, typename glue_type>
inline
void
op_scalar_div_post::apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_type>, op_scalar_div_post>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  out = in.m;  // implicit conversion to 'Mat<eT>'
  
        eT* out_mem = out.memptr();
  const eT  k       = in.aux;
  const u32 n_elem  = out.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] /= k;
    }
  }



template<typename T1, typename T2>
inline
void
op_scalar_div_post::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_plus>, op_scalar_div_post>& in)
  {
  arma_extra_debug_sigprint();
  
  unwrap<T1> tmp1(in.m.A);
  unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix addition");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A_mem[i] + B_mem[i]) / k;
    }
  
  }



template<typename T1, typename T2>
inline
void
op_scalar_div_post::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_minus>, op_scalar_div_post>& in)
  {
  arma_extra_debug_sigprint();
  
  unwrap<T1> tmp1(in.m.A);
  unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "matrix subtraction");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A_mem[i] - B_mem[i]) / k;
    }
  
  }



template<typename T1, typename T2>
inline
void
op_scalar_div_post::apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2,glue_schur>, op_scalar_div_post>& in)
  {
  arma_extra_debug_sigprint();
  
  unwrap<T1> tmp1(in.m.A);
  unwrap<T2> tmp2(in.m.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_assert_same_size(A, B, "schur product");
  
  // no alias problems
  out.set_size(A.n_rows, A.n_cols);
  
        eT* out_mem = out.memptr();
  const eT* A_mem   = A.mem;
  const eT* B_mem   = B.mem;
  
  const eT  k       = in.aux;
  const u32 n_elem  = A.n_elem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_mem[i] = (A_mem[i] * B_mem[i]) / k;
    }
  
  }



//! @}
