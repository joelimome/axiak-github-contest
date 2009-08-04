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


//! \addtogroup fn_norm
//! @{

template<typename T1>
inline
typename T1::elem_type
norm(const Base<typename T1::elem_type,T1>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> A_tmp(X.get_ref());
  const Mat<eT>& A = A_tmp.M;

  arma_debug_check(    (A.n_elem == 0),                      "norm(): given object has no elements" );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): first argument must be a vector" );
  arma_debug_check(    (k == 0),                             "norm(): second argument must be greater than zero" );


  if(k==1)
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      acc += std::abs(A.mem[i]);
      }
    
    return acc;
    }
  else
  if(k==2)
    {
    if(is_complex<eT>::value == false)
      {
      eT acc = eT(0);
      
      for(u32 i=0; i<A.n_elem; ++i)
        {
        const eT tmp = A.mem[i];
        acc += tmp*tmp;
        }
      
      return std::sqrt(acc);
      }
    else
      {
      eT acc = eT(0);
      
      for(u32 i=0; i<A.n_elem; ++i)
        {
        acc += std::abs(A.mem[i]);
        }
      
      return std::sqrt(acc);
      
      }
    }
  else
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      acc += std::pow(std::abs(A.mem[i]), int(k));
      }
    
    return std::pow(acc, eT(1)/eT(k));
    }
  
  }



//
// giving vector arguments works, as operator-() takes two mat arguments.
// i.e. it has no specific form for dealing with colvec and rowvec,
// and colvec and rowvec are derived from mat (hence they are a type of mat)

template<typename T1, typename T2>
inline
typename T1::elem_type
norm(const Glue<T1,T2,glue_minus>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type,typename T2::elem_type>::check();
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check(  ( (A.n_elem == 0) || (B.n_elem == 0) ), "norm(): one or more of given objects has no elements" );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): non-vector argument detected" );
  arma_debug_check( !( (B.n_rows == 1) || (B.n_cols == 1) ), "norm(): non-vector argument detected" );
  arma_debug_check(    (A.n_elem != B.n_elem),               "norm(): vectors have different lengths" );
  arma_debug_check(    (k == 0),                             "norm(): second argument must be greater than zero" );
  
  if(k==1)
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      acc += std::abs(A.mem[i] - B.mem[i]);
      }
    
    return acc;
    }
  else
  if(k==2)
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      const eT tmp = A.mem[i] - B.mem[i];
      
      acc += tmp*tmp;
      }
    
    return std::sqrt(acc);
    }
  else
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      acc += std::pow( std::abs(A.mem[i] - B.mem[i]), int(k) );
      }
    
    return std::pow(acc, eT(1)/eT(k));
    }
  
  }



template<typename T1, typename T2>
inline
typename T1::elem_type
norm(const Glue<T1,T2,glue_plus>& X, const u32 k)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<typename T1::elem_type,typename T2::elem_type>::check();
  
  const unwrap<T1> tmp1(X.A);
  const unwrap<T2> tmp2(X.B);
  
  typedef typename T1::elem_type eT;
  
  const Mat<eT>& A = tmp1.M;
  const Mat<eT>& B = tmp2.M;
  
  arma_debug_check(  ( (A.n_elem == 0) || (B.n_elem == 0) ), "norm(): one or more of given objects has no elements" );
  arma_debug_check( !( (A.n_rows == 1) || (A.n_cols == 1) ), "norm(): non-vector argument detected" );
  arma_debug_check( !( (B.n_rows == 1) || (B.n_cols == 1) ), "norm(): non-vector argument detected" );
  arma_debug_check(    (A.n_elem != B.n_elem),               "norm(): vectors have different lengths" );
  arma_debug_check(    (k == 0),                             "norm(): second argument must be greater than zero" );
  
  if(k==1)
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      acc += std::abs(A.mem[i] + B.mem[i]);
      }
    
    return acc;
    }
  else
  if(k==2)
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      const eT tmp = A.mem[i] + B.mem[i];
      
      acc += tmp*tmp;
      }
    
    return std::sqrt(acc);
    }
  else
    {
    eT acc = eT(0);
    
    for(u32 i=0; i<A.n_elem; ++i)
      {
      acc += std::pow( std::abs(A.mem[i] + B.mem[i]), int(k) );
      }
    
    return std::pow(acc, eT(1)/eT(k));
    }
  
  }


//! @}
