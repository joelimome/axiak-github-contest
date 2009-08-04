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


//! \addtogroup operator_relational
//! @{



template<typename eT, typename T1, typename T2>
inline
umat
operator==
(const Base<eT,T1>& X, const Base<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
    
  arma_debug_assert_same_size(A, B, "operator==");
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  const eT* B_mem = B.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] == B_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename eT, typename T1, typename T2>
inline
umat
operator!=
(const Base<eT,T1>& X, const Base<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
    
  arma_debug_assert_same_size(A, B, "operator!=");
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  const eT* B_mem = B.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] != B_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename eT, typename T1, typename T2>
inline
umat
operator>=
(const Base<eT,T1>& X, const Base<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
    
  arma_debug_assert_same_size(A, B, "operator>=");
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  const eT* B_mem = B.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] >= B_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename eT, typename T1, typename T2>
inline
umat
operator<=
(const Base<eT,T1>& X, const Base<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
    
  arma_debug_assert_same_size(A, B, "operator<=");
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  const eT* B_mem = B.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] == B_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename eT, typename T1, typename T2>
inline
umat
operator>
(const Base<eT,T1>& X, const Base<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
    
  arma_debug_assert_same_size(A, B, "operator>");
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  const eT* B_mem = B.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] > B_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename eT, typename T1, typename T2>
inline
umat
operator<
(const Base<eT,T1>& X, const Base<eT,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
    
  arma_debug_assert_same_size(A, B, "operator<");
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  const eT* B_mem = B.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] < B_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator==
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] == val)
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator!=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] != val)
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator>=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] >= val)
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<=
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] <= val)
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator>
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] > val)
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type val)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(A_mem[i] < val)
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator==
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  return operator==(X,val);
  }



template<typename T1>
inline
umat
operator!=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  return operator!=(X,val);
  }



template<typename T1>
inline
umat
operator>=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val >= A_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<=
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val <= A_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator>
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
    
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val > A_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



template<typename T1>
inline
umat
operator<
(const typename T1::elem_type val, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp1(X.get_ref());
  const Mat<eT>& A = tmp1.M;
  
  umat out(A.n_rows, A.n_cols);
  
  const eT* A_mem = A.mem;
  
  typedef typename umat::elem_type umat_elem_type;
  umat_elem_type* out_mem = out.memptr();
  
  for(u32 i=0; i<A.n_elem; ++i)
    {
    if(val < A_mem[i])
      {
      out_mem[i] = umat_elem_type(1);
      }
    else
      {
      out_mem[i] = umat_elem_type(0);
      }
    }
  
  return out;
  }



//! @}
