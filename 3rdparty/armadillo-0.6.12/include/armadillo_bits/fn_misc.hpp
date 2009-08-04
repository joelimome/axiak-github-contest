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


//! \addtogroup fn_misc
//! @{

//! \brief
//! Generate a vector with 'num' elements.
//! The values of the elements linearly increase from 'start' upto (and including) 'end'.

template<typename eT>
inline
Mat<eT>
linspace(const eT start, const eT end, const u32 num, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (num < 2), "linspace(): num must be >= 2");
  
  Mat<eT> x;
  
  if(dim == 0)
    {
    x.set_size(num,1);  // column vector
    }
  else
    {
    x.set_size(1,num);  // row vector
    }
  
  
  const eT delta = (end-start)/(num-1);
  
  x[0] = start;
  
  for(u32 i=1; i<num; ++i)
    {
    x[i] = x[i-1] + delta;
    }
  
  return x; 
  }



inline
mat
linspace(const double start, const double end, const u32 num, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();
  return linspace<double>(start, end, num, dim);
  }






//
// reshape

template<typename T1>
inline
Mat<typename T1::elem_type>
reshape(const Base<typename T1::elem_type,T1>& X, const u32 in_n_rows, const u32 in_n_cols, const u32 dim = 0)
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  Mat<eT> out;

  const unwrap<T1> A_tmp(X.get_ref());
  const Mat<eT>& A = A_tmp.M;

  const u32 in_n_elem = in_n_rows * in_n_cols;

  arma_debug_check( (A.n_elem != in_n_elem), "reshape(): incompatible dimensions");
  arma_debug_check( (dim > 1), "reshape(): dim must be 0 or 1");

  if(dim == 0)
    {
    out = A;

    access::rw(out.n_rows) = in_n_rows;
    access::rw(out.n_cols) = in_n_cols;

    return out;
    }
  else
    {
    out.set_size(in_n_rows, in_n_cols);
    
    eT* out_mem = out.memptr();
    u32 i = 0;
    
    for(u32 row=0; row<A.n_rows; ++row)
      {
      for(u32 col=0; col<A.n_cols; ++col)
        {
        out_mem[i] = A.at(row,col);
        ++i;
        }
      }
    
    return out;
    }
  }



//
// real

template<typename T, typename T1>
inline
Mat<T>
real(const Base<std::complex<T>, T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;

  const unwrap<T1> A_tmp(X.get_ref());
  const Mat<eT>& A = A_tmp.M;
  
  Mat<T> out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  T* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::real(A_mem[i]);
    }
  
  return out;
  }



//
// imag

template<typename T, typename T1>
inline
Mat<T>
imag(const Base<std::complex<T>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef std::complex<T> eT;

  const unwrap<T1> A_tmp(X.get_ref());
  const Mat<eT>& A = A_tmp.M;
  
  Mat<T> out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  T* out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::imag(A_mem[i]);
    }
  
  return out;
  }



//
// log_add

template<typename eT>
inline
eT
log_add(eT log_a, eT log_b)
  {
  if(log_a < log_b)
    {
    std::swap(log_a, log_b);
    }
  
  const eT negdelta = log_b - log_a;
  
  if( (negdelta < Math<eT>::log_min()) || (arma_isfinite(negdelta) == false) )
    {
    return log_a;
    }
  else
    {
    #if defined(ARMA_HAVE_LOG1P)
      return (log_a + log1p(std::exp(negdelta)));
    #else
      return (log_a + std::log(1.0 + std::exp(negdelta)));
    #endif
    }
  }



template<typename eT>
inline
eT
trunc_log(const eT x)
  {
  if(std::numeric_limits<eT>::is_iec559)
    {
    if(x == std::numeric_limits<eT>::infinity())
      {
      return Math<eT>::log_max();
      }
    if(x <= 0)
      {
      return Math<eT>::log_min();
      }
    }
  else
    {
    return std::log(x);
    }
  }



template<typename eT>
inline
eT
trunc_exp(const eT x)
  {
  if(std::numeric_limits<eT>::is_iec559 && (x >= Math<eT>::log_max() ))
    {
    return std::numeric_limits<eT>::max();
    }
  else
    {
    return std::exp(x);
    }
  }



//
// log

template<typename T1>
inline
const Op<T1, op_log>
log(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_log>(A.get_ref());
  }



//
// trunc_log

template<typename T1>
inline
const Op<T1, op_trunc_log>
trunc_log(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_trunc_log>(A.get_ref());
  }



//
// log10

template<typename T1>
inline
const Op<T1, op_log10>
log10(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_log10>(A.get_ref());
  }



//
// exp

template<typename T1>
inline
const Op<T1, op_exp>
exp(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_exp>(A.get_ref());
  }



//
// trunc_exp

template<typename T1>
inline
const Op<T1, op_trunc_exp>
trunc_exp(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_trunc_exp>(A.get_ref());
  }



//
// abs

template<typename T1>
inline
Mat<typename T1::pod_type>
abs(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> A_tmp(X.get_ref());

  // if T1 is a complex matrix,
  // pod_type is the underlying type used by std::complex;
  // otherwise pod_type is the same as elem_type
  
  typedef typename T1::elem_type  in_eT;
  typedef typename T1::pod_type  out_eT;

  const Mat<in_eT>& A = A_tmp.M;
  
  Mat<out_eT> out(A.n_rows, A.n_cols);
  
  const in_eT* A_mem   = A.mem;
  out_eT*      out_mem = out.memptr();
  
  for(u32 i=0; i<out.n_elem; ++i)
    {
    out_mem[i] = std::abs(A_mem[i]);
    }
  
  return out;
  }



//
// fabs

template<typename T1>
inline
Mat<typename T1::pod_type>
fabs(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return abs(A);
  }



//
// square

template<typename T1>
inline
const Op<T1, op_square>
square(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_square>(A.get_ref());
  }



//
// sqrt

template<typename T1>
inline
const Op<T1, op_sqrt>
sqrt(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sqrt>(A.get_ref());
  }



// pow

template<typename T1>
inline
const Op<T1, op_pow>
pow(const Base<typename T1::elem_type,T1>& A, const typename T1::elem_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_pow>(A.get_ref(), exponent);
  }



// pow, specialised handling (non-complex exponent for complex matrices)

template<typename T1>
inline
const Op<T1, op_pow>
pow(const Base<typename T1::elem_type,T1>& A, const typename T1::elem_type::value_type exponent)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_pow>(A.get_ref(), eT(exponent));
  }



// pow_s32  (integer exponent)

template<typename T1>
inline
const Op<T1, op_pow_s32>
pow(const Base<typename T1::elem_type,T1>& A, const s32 exponent)
  {
  arma_extra_debug_sigprint();
  
  if(exponent >= 0)
    {
    return Op<T1, op_pow_s32>(A.get_ref(), exponent, 0);
    }
  else
    {
    return Op<T1, op_pow_s32>(A.get_ref(), -exponent, 1);
    }
  }



// conj

template<typename T, typename T1>
inline
const Op<T1, op_conj>
conj(const Base<std::complex<T>,T1>& A)
  {
  arma_extra_debug_sigprint();

  return Op<T1, op_conj>(A.get_ref());
  }



template<typename T1>
inline
const T1&
conj(const Op<T1, op_conj>& A)
  {
  arma_extra_debug_sigprint();
  
  return A.m;
  }



//! the conjugate of the transpose of a complex matrix is the same as the hermitian transpose
template<typename T1>
inline
const Op<T1, op_htrans>
conj(const Op<T1, op_trans>& A)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< is_complex<typename T1::elem_type>::value == false >::apply();

  return Op<T1, op_htrans>(A.m);
  }


//! @}
