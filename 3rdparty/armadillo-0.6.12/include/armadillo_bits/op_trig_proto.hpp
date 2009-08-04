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

class op_cos
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cos>& in);
  };



class op_acos
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_acos>& in);
  template<typename T, typename T1> inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_acos>& in);
  };



class op_cosh
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cosh>& in);
  };
  


class op_acosh
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_acosh>& in);
  template<typename T, typename T1> inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_acosh>& in);
  };
  


// sin family

class op_sin
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sin>& in);
  };



class op_asin
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_asin>& in);
  template<typename T, typename T1> inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_asin>& in);
  };



class op_sinh
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_sinh>& in);
  };
  


class op_asinh
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_asinh>& in);
  };



// tan family

class op_tan
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_tan>& in);
  };



class op_atan
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_atan>& in);
  template<typename T, typename T1> inline static void apply(Mat< std::complex<T> >& out, const Op<T1,op_atan>& in);
  };



class op_tanh
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_tanh>& in);
  };
  


class op_atanh
  {
  public:
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_atanh>& in);
  };


//! @}
