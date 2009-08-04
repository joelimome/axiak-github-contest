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



//! \addtogroup op_trig
//! @{


//
// trigonometric functions:
// cos family: cos, acos, cosh, acosh
// sin family: sin, asin, sinh, asinh
// tan family: tan, atan, tanh, atanh

// cos family

template<typename T1>
inline
void
op_cos::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cos>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::cos(A.mem[i]);
    }
  
  }



template<typename T1>
inline
void
op_acos::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_acos>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::acos(A.mem[i]);
    }
  
  }



template<typename T, typename T1>
inline
void
op_acos::apply(Mat< std::complex<T> >& out, const Op<T1,op_acos>& in)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST)
    {
    typedef typename std::complex<T> eT;
    isnt_same_type<eT, typename T1::elem_type>::check();
    
    const unwrap<T1> tmp(in.m);
    
    const Mat<eT>& A = tmp.M;
    const u32 n_elem = A.n_elem;
    
    out.set_size(A.n_rows, A.n_cols);
    eT* out_ptr = out.memptr();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = boost::math::acos(A.mem[i]);
      }
    }
  #else
    {
    arma_stop("op_acos::apply(): need Boost libraries");
    }
  #endif
  }



template<typename T1>
inline
void
op_cosh::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cosh>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::cosh(A.mem[i]);
    }
  
  }



template<typename T1>
inline
void
op_acosh::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_acosh>& in)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST)
    {
    typedef typename T1::elem_type eT;
    
    const unwrap<T1> tmp(in.m);
    
    const Mat<eT>& A = tmp.M;
    const u32 n_elem = A.n_elem;
    
    out.set_size(A.n_rows, A.n_cols);
    eT* out_ptr = out.memptr();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = boost::math::acosh(A.mem[i]);
      }
    }
  #else
    {
    arma_stop("op_acosh::apply(): need Boost libraries");
    }
  #endif
  
  }



// sin family

template<typename T1>
inline
void
op_sin::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sin>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::sin(A.mem[i]);
    }
  
  }



template<typename T1>
inline
void
op_asin::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_asin>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::acos(A.mem[i]);
    }
  
  }



template<typename T, typename T1>
inline
void
op_asin::apply(Mat< std::complex<T> >& out, const Op<T1,op_asin>& in)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST)
    {
    typedef typename std::complex<T> eT;
    isnt_same_type<eT, typename T1::elem_type>::check();
    
    const unwrap<T1> tmp(in.m);
    
    const Mat<eT>& A = tmp.M;
    const u32 n_elem = A.n_elem;
    
    out.set_size(A.n_rows, A.n_cols);
    eT* out_ptr = out.memptr();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = boost::math::asin(A.mem[i]);
      }
    }
  #else
    {
    arma_stop("op_asin::apply(): need Boost libraries");
    }
  #endif
  }



template<typename T1>
inline
void
op_sinh::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sinh>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::sinh(A.mem[i]);
    }
  
  }



template<typename T1>
inline
void
op_asinh::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_asinh>& in)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST)
    {
    typedef typename T1::elem_type eT;
    
    const unwrap<T1> tmp(in.m);
  
    const Mat<eT>& A = tmp.M;
    const u32 n_elem = A.n_elem;
    
    out.set_size(A.n_rows, A.n_cols);
    eT* out_ptr = out.memptr();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = boost::math::asinh(A.mem[i]);
      }
    }
  #else
    {
    arma_stop("op_asinh::apply(): need Boost libraries");
    }
  #endif
  }



// tan family

template<typename T1>
inline
void
op_tan::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_tan>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::tan(A.mem[i]);
    }
  
  }



template<typename T1>
inline
void
op_atan::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_atan>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::atan(A.mem[i]);
    }
  
  }



template<typename T, typename T1>
inline
void
op_atan::apply(Mat< std::complex<T> >& out, const Op<T1,op_atan>& in)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST)
    {
    typedef typename std::complex<T> eT;
    isnt_same_type<eT, typename T1::elem_type>::check();
    
    const unwrap<T1> tmp(in.m);
    
    const Mat<eT>& A = tmp.M;
    const u32 n_elem = A.n_elem;
    
    out.set_size(A.n_rows, A.n_cols);
    eT* out_ptr = out.memptr();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = boost::math::atan(A.mem[i]);
      }
    }
  #else
    {
    arma_stop("op_asin::apply(): need Boost libraries");
    }
  #endif
  }



template<typename T1>
inline
void
op_tanh::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_tanh>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(in.m);
  
  const Mat<eT>& A = tmp.M;
  const u32 n_elem = A.n_elem;
  
  out.set_size(A.n_rows, A.n_cols);
  eT* out_ptr = out.memptr();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    out_ptr[i] = std::tanh(A.mem[i]);
    }
  
  }



template<typename T1>
inline
void
op_atanh::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_atanh>& in)
  {
  arma_extra_debug_sigprint();
  
  #if defined(ARMA_USE_BOOST)
    {
    typedef typename T1::elem_type eT;
    
    const unwrap<T1> tmp(in.m);
  
    const Mat<eT>& A = tmp.M;
    const u32 n_elem = A.n_elem;
    
    out.set_size(A.n_rows, A.n_cols);
    eT* out_ptr = out.memptr();
    
    for(u32 i=0; i<n_elem; ++i)
      {
      out_ptr[i] = boost::math::atanh(A.mem[i]);
      }
    }
  #else
    {
    arma_stop("op_atanh::apply(): need Boost libraries");
    }
  #endif
  
  }



//! @}
