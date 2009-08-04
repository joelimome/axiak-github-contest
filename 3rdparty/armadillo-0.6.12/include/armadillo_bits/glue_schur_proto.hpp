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


//! \addtogroup glue_schur
//! @{


//! Class which implements the immediate Schur product (element-wise multiplication) of two or more matrices
class glue_schur
  {
  public:
  
  // mat
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue<Mat<eT>,Mat<eT>,glue_schur>& X);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue< Glue<Mat<eT>,Mat<eT>,glue_schur>, Mat<eT>, glue_schur>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_schur>& X);

  
  // mat, inplace
  
  template<typename eT> inline static void apply_inplace(Mat<eT>& out, const Mat<eT>& B);
  
  template<typename T1, typename op_type>
  inline static void apply_inplace(Mat<typename T1::elem_type>& out, const Op<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply_inplace(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_type>& X);

  
  // matrix schur product with different element types
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(Mat<typename promote_type<eT1,eT2>::result>& out, const Mat<eT1>& X, const Mat<eT2>& Y);
  
  };



class glue_schur_diag
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const T1& A, const Op<T2,op_diagmat>& B);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& A, const Op<T2,op_diagmat>& B);
  
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1, Op<T2,op_diagmat>, glue_schur_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, T2, glue_schur_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_schur_diag>& X);
  
  };

//! @}

