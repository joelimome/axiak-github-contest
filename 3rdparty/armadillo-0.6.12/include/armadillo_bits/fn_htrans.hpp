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


//! \addtogroup fn_htrans
//! @{


template<typename T>
inline
const Op<Mat< std::complex<T> >, op_htrans>
htrans(const Mat< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Mat< std::complex<T> >, op_htrans>(X);
  }



template<typename T>
inline
const Op<Row< std::complex<T> >, op_htrans>
htrans(const Row< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Row< std::complex<T> >, op_htrans>(X);
  }



template<typename T>
inline
const Op<Col< std::complex<T> >, op_htrans>
htrans(const Col< std::complex<T> >& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<Col< std::complex<T> >, op_htrans>(X);
  }



template<typename T, typename T1>
inline
const Op<T1, op_htrans>
htrans(const Base<std::complex<T>,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_htrans>(X.get_ref());
  }



//! two consecutive hermitian transpose operations cancel each other
template<typename T1>
inline
const T1&
htrans(const Op<T1, op_htrans>& X)
  {
  arma_extra_debug_sigprint();
  arma_extra_debug_print("htrans(): removing op_htrans");
  
  return X.m;
  }


//! @}
