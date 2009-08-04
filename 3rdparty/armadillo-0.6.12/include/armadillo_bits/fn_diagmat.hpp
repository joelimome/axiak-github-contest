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


//! \addtogroup fn_diagmat
//! @{


//! interpret a mat as a diagonal matrix (i.e. off-diagonal entries are zero)
template<typename eT, typename T1>
inline
const Op<T1, op_diagmat>
diagmat(const Base<eT,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_diagmat>(X.get_ref());
  }



//! interpret a colvec as a diagonal matrix
template<typename eT>
inline
const Op<Mat<eT>, op_diagmat_vec>
diagmat(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Mat<eT>, op_diagmat_vec>(X);
  }



//! interpret a rowvec as a diagonal matrix
template<typename eT>
inline
const Op<Mat<eT>, op_diagmat_vec>
diagmat(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Mat<eT>, op_diagmat_vec>(X);
  }



//! create a diagonal matrix out of subview_col
template<typename eT>
inline
Mat<eT>
diagmat(const subview_col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.n_elem, X.n_elem);
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out.at(i,i) = X[i];
    }
  
  return out;
  }



//! create a diagonal matrix out of subview_row
template<typename eT>
inline
Mat<eT>
diagmat(const subview_row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.n_elem, X.n_elem);
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out.at(i,i) = X[i];
    }
  
  return out;
  }



//! create a diagonal matrix out of diagview
template<typename eT>
inline
Mat<eT>
diagmat(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.n_elem, X.n_elem);
  
  for(u32 i=0; i<X.n_elem; ++i)
    {
    out.at(i,i) = X[i];
    }
  
  return out;
  }



//! two consecutive diagmat operations are equivalent to one diagmat operation
template<typename T1>
inline
const Op<T1, op_diagmat>&
diagmat(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("diagmat(): two consecutive diagmat() operations detected");
  
  return X;
  }



//! @}
