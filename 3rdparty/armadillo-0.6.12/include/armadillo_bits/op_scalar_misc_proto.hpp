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



//! \addtogroup op_scalar_misc
//! @{


//! 'add scalar to matrix' operation
class op_scalar_plus
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_plus>& in);
  };



//! 'subtract matrix from a scalar' operation
class op_scalar_minus_pre
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_minus_pre>& in);
  };



//! 'subtract scalar from matrix' operation
class op_scalar_minus_post
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_minus_post>& in);
  };


//! 'multiply matrix by a scalar' operation
class op_scalar_times
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_times>& in);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_type>,  op_scalar_times>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_plus>,  op_scalar_times>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_minus>, op_scalar_times>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_schur>, op_scalar_times>& in);
  
  
  };



//! 'divide scalar by a matrix' operation
class op_scalar_div_pre
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_div_pre>& in);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Op<Mat<eT>,op_scalar_div_pre>& in);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_type>, op_scalar_div_pre>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_plus>,  op_scalar_div_pre>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_minus>, op_scalar_div_pre>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_schur>, op_scalar_div_pre>& in);
  
  };



//! 'divide matrix by a scalar' operation
class op_scalar_div_post
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_scalar_div_post>& in);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Op<Mat<eT>,op_scalar_div_post>& in);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_type>, op_scalar_div_post>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_plus>,  op_scalar_div_post>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_minus>, op_scalar_div_post>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<Glue<T1,T2,glue_schur>, op_scalar_div_post>& in);
  
  };



//! @}
