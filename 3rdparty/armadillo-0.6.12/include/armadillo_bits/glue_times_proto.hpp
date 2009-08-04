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


//! \addtogroup glue_times
//! @{


//! Class which implements the immediate multiplication of two or more matrices
class glue_times
  {
  public:
  
  
  template<typename eT>
  arma_inline static u32  mul_storage_cost(const Mat<eT>& X, const Mat<eT>& Y);
  
  template<typename eT>
  inline static void apply_noalias(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const Mat<eT>& C);
  
  

  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X);

  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue<Mat<eT>,Mat<eT>,glue_times>& X);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue< Glue<Mat<eT>,Mat<eT>, glue_times>, Mat<eT>, glue_times>& X);
  
  
  
  template<typename eT>
  inline static void apply_inplace(Mat<eT>& out, const Mat<eT>& B);
  
  template<typename T1, typename op_type>
  inline static void apply_inplace(Mat<typename T1::elem_type>& out, const Op<T1, op_type>& X);
  
  template<typename T1, typename T2, typename glue_type>
  inline static void apply_inplace(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_type>& X);
  
  
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,Op<T2,op_trans>,glue_times >& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_trans>,T2,glue_times>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_trans>,Op<T2,op_trans>,glue_times>& X);
  
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue< Op<T1, op_neg>, T2, glue_times>& X);
  
  
  template<typename eT>
  inline static eT direct_rowvec_mat_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);

  template<typename eT>
  inline static eT direct_rowvec_diagmat_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagmat_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  template<typename eT>
  inline static eT direct_rowvec_invdiagvec_colvec(const eT* A_mem, const Mat<eT>& B, const eT* C_mem);
  
  
  // matrix multiplication with different element types
  
  template<typename eT1, typename eT2>
  inline static void apply_mixed(Mat<typename promote_type<eT1,eT2>::result>& out, const Mat<eT1>& X, const Mat<eT2>& Y);
  
  };



class glue_times_diag
  {
  public:
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const T1& A, const Op<T2,op_diagmat>& B);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& A, const T2& B);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& A, const Op<T2,op_diagmat>& B);
  
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1, Op<T2,op_diagmat>, glue_times_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, T2, glue_times_diag>& X);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_times_diag>& X);
  
  };



class glue_times_vec
  {
  public:
  
  template<typename eT>
  inline static void mul_col_row(Mat<eT>& out, const eT* A_mem, const eT* B_mem);

  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1, Col<typename T1::elem_type>, glue_times_vec>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1, Row<typename T1::elem_type>, glue_times_vec>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Col<typename T1::elem_type>, T1, glue_times_vec>& X);
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Row<typename T1::elem_type>, T1, glue_times_vec>& X);
  
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue<Col<eT>, Row<eT>, glue_times_vec>& X);

  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue< Op<Row<eT>, op_trans>, Row<eT>, glue_times_vec>& X);
  
  template<typename eT>
  inline static void apply(Mat<eT>& out, const Glue< Col<eT>, Op<Col<eT>, op_trans>, glue_times_vec>& X);
  
  
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Glue<Op<T1, op_trans>, Col<typename T1::elem_type>, glue_times_vec>& X);
  
  };


//! @}

