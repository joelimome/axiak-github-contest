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


//! \addtogroup op_diagmat
//! @{


//! convert a mat/rowvec/colvec to a diagonal matrix

class op_diagmat
  {
  private:

  template<typename eT>
  inline static void zero_offdiag(Mat<eT>& X);
  
  
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_div>,   op_diagmat>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_minus>, op_diagmat>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_plus>,  op_diagmat>& in);
  
  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_schur>, op_diagmat>& in);

  template<typename T1, typename T2>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op< Glue<T1,T2, glue_times>, op_diagmat>& in);
  
  };


class op_diagmat_vec
  {
  public:
  
  template<typename T1>
  inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_diagmat_vec>& in);
  };


//! @}
