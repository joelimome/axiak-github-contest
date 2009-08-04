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


//! \addtogroup op_dot
//! @{



//! for two arrays
template<typename eT>
inline
arma_pure
eT
op_dot::direct_dot(const u32 n_elem, const eT* const A, const eT* const B)
  {
  arma_extra_debug_sigprint();
  
  eT val1 = eT(0);
  eT val2 = eT(0);
  
  u32 i,j;
  for(i=0, j=1; j<n_elem; i+=2, j+=2)
    {
    val1 += A[i] * B[i];
    val2 += A[j] * B[j];
    }
  
  if(i < n_elem)
    {
    val1 += A[i] * B[i];
    }
  
  return val1+val2;
  }



//! for three arrays
template<typename eT>
inline
arma_pure
eT
op_dot::direct_dot(const u32 n_elem, const eT* const A, const eT* const B, const eT* C)
  {
  arma_extra_debug_sigprint();
  
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += A[i] * B[i] * C[i];
    }

  return val;
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
op_dot::apply(const Base<typename T1::elem_type,T1>& A_orig, const Base<typename T1::elem_type,T2>& B_orig)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_to_elem_access<T1> A(A_orig.get_ref());
  const unwrap_to_elem_access<T2> B(B_orig.get_ref());

  arma_debug_check( (A.M.n_elem != B.M.n_elem), "dot(): objects must have the same number of elements" );
  
  const u32 n_elem = A.M.n_elem;
  eT val = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    val += A[i] * B[i];
    }
  
  return val;
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
op_norm_dot::apply(const Base<typename T1::elem_type,T1>& A_orig, const Base<typename T1::elem_type,T2>& B_orig)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_to_elem_access<T1> A(A_orig.get_ref());
  const unwrap_to_elem_access<T2> B(B_orig.get_ref());

  arma_debug_check( (A.M.n_elem != B.M.n_elem), "norm_dot(): objects must have the same number of elements" );
  
  const u32 n_elem = A.M.n_elem;
  
  eT acc1 = eT(0);
  eT acc2 = eT(0);
  eT acc3 = eT(0);
  
  for(u32 i=0; i<n_elem; ++i)
    {
    const eT tmpA = A[i];
    const eT tmpB = B[i];
    
    acc1 += tmpA * tmpA;
    acc2 += tmpB * tmpB;
    acc3 += tmpA * tmpB;
    }
    
  return acc3 / ( std::sqrt(acc1 * acc2) );   // TODO: this only makes sense for eT = float, double or complex
  }


//! @}
