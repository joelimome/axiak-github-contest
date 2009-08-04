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


//! \addtogroup fn_trig
//! @{

//
// trigonometric functions:
// cos family: cos, acos, cosh, acosh
// sin family: sin, asin, sinh, asinh
// tan family: tan, atan, tanh, atanh


//
// cos

template<typename T1>
inline
const Op<T1, op_cos>
cos(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_cos>(A.get_ref());
  }



//
// acos

template<typename T1>
inline
const Op<T1, op_acos>
acos(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_acos>(A.get_ref());
  }



//
// cosh

template<typename T1>
inline
const Op<T1, op_cosh>
cosh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_cosh>(A.get_ref());
  }



//
// acosh

template<typename T1>
inline
const Op<T1, op_acosh>
acosh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_acosh>(A.get_ref());
  }



//
// sin

template<typename T1>
inline
const Op<T1, op_sin>
sin(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sin>(A.get_ref());
  }



//
// asin

template<typename T1>
inline
const Op<T1, op_asin>
asin(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_asin>(A.get_ref());
  }



//
// sinh

template<typename T1>
inline
const Op<T1, op_sinh>
sinh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_sinh>(A.get_ref());
  }



//
// asinh

template<typename T1>
inline
const Op<T1, op_asinh>
asinh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_asinh>(A.get_ref());
  }



//
// tan

template<typename T1>
inline
const Op<T1, op_tan>
tan(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_tan>(A.get_ref());
  }



//
// atan

template<typename T1>
inline
const Op<T1, op_atan>
atan(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_atan>(A.get_ref());
  }



//
// tanh

template<typename T1>
inline
const Op<T1, op_tanh>
tanh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_tanh>(A.get_ref());
  }



//
// atanh

template<typename T1>
inline
const Op<T1, op_atanh>
atanh(const Base<typename T1::elem_type,T1>& A)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_atanh>(A.get_ref());
  }



//! @}
