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


//! \addtogroup Row
//! @{

//! Class for row vectors (matrices with only one row)

template<typename eT>
class Row : public Mat<eT>, public Base_vec<eT, Row<eT> >
  {
  public:
  
  typedef eT elem_type;
  typedef typename get_pod_type<elem_type>::pod_type pod_type;
  
  
  inline                     Row();
  inline explicit            Row(const u32 N);
  
  inline                     Row(const char* text);
  inline const Row&    operator=(const char* text);  // TODO: std::string input
  
  inline                     Row(const Row& X);
  inline const Row&    operator=(const Row& X);
  inline const Row&   operator*=(const Row& X);
  
  //inline explicit            Row(const Mat<eT>& X);
  inline                     Row(const Mat<eT>& X);
  inline const Row&    operator=(const Mat<eT>& X);
  inline const Row&   operator*=(const Mat<eT>& X);
  
  inline                     Row(const eT* aux_mem, const u32 aux_length);

  template<typename T1, typename T2>
  inline explicit Row(const Base<pod_type,T1>& A, const Base<pod_type,T2>& B);

  inline                     Row(const subview<eT>& X);
  inline const Row&    operator=(const subview<eT>& X);
  inline const Row&   operator*=(const subview<eT>& X);
  
  inline explicit            Row(const diagview<eT>& X);
  inline const Row&    operator=(const diagview<eT>& X);
  inline const Row&   operator*=(const diagview<eT>& X);
  
  template<typename T1, typename op_type> inline                   Row(const Op<T1, op_type> &X);
  template<typename T1, typename op_type> inline const Row&  operator=(const Op<T1, op_type> &X);
  template<typename T1, typename op_type> inline const Row& operator*=(const Op<T1, op_type> &X);
  
  template<typename T1, typename T2, typename glue_type> inline                   Row(const Glue<T1, T2, glue_type> &X);
  template<typename T1, typename T2, typename glue_type> inline const Row&  operator=(const Glue<T1, T2, glue_type> &X);
  template<typename T1, typename T2, typename glue_type> inline const Row& operator*=(const Glue<T1, T2, glue_type> &X);
  
  
  inline void set_size(const u32 N);
  inline void set_size(const u32 n_rows, const u32 n_cols);
  
  inline void zeros();
  inline void zeros(const u32 N);
  inline void zeros(const u32 n_rows, const u32 n_cols);


  inline void load(const std::string name, const file_type type = auto_detect);
  };


//! @}
