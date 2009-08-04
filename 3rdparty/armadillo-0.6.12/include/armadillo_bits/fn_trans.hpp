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


//! \addtogroup fn_trans
//! @{


template<typename T1>
inline
const Op<T1, op_trans>
trans(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_trans>(X.get_ref());
  }



template<typename eT>
inline
const Op<Row<eT>, op_trans>
trans(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Row<eT>, op_trans>(X);
  }



template<typename eT>
inline
const Op<Col<eT>, op_trans>
trans(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Col<eT>, op_trans>(X);
  }



//! two consecutive transpose operations cancel each other
template<typename T1>
inline
const T1&
trans(const Op<T1, op_trans>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("trans(): removing op_trans");
  
  return X.m;
  }



//! transpose of a diagonal matrix gives you the original matrix back
template<typename T1>
inline
const Op<T1, op_diagmat>
trans(const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("trans(): not applying op_trans to diagonal matrix");
  
  return X;
  }



//! the transpose of the conjugate of a complex matrix is equivalent to the hermitian transpose
template<typename T1>
inline
const Op<T1, op_htrans>
trans(const Op<T1, op_conj>& X)
  {
  arma_extra_debug_sigprint();

  arma_type_check< is_complex<typename T1::elem_type>::value == false >::apply();

  return Op<T1, op_htrans>(X.m);
  }



//! @}
