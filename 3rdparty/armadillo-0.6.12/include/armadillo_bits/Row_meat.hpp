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

template<typename eT>
inline
Row<eT>::Row()
  : Mat<eT>()
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
Row<eT>::Row(const u32 in_n_elem)
  : Mat<eT>(1,in_n_elem)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
Row<eT>::Row(const char* text)
  : Mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }
  


template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
Row<eT>::Row(const Row<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const Row<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  return *this;
  }



template<typename eT>
inline Row<eT>::Row(const Mat<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



//! construct a row vector from a given auxillary array
template<typename eT>
inline
Row<eT>::Row(const eT* aux_mem, const u32 aux_length)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size(1, aux_length);
  arma_check( (Mat<eT>::n_elem != aux_length), "Row(): don't know how to handle the given array" );

  syslib::copy_elem( Mat<eT>::memptr(), aux_mem, Mat<eT>::n_elem );
  }



template<typename eT>
template<typename T1, typename T2>
inline
Row<eT>::Row
  (
  const Base<typename Row<eT>::pod_type, T1>& A,
  const Base<typename Row<eT>::pod_type, T2>& B
  )
  : Mat<eT>(A,B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
Row<eT>::Row(const subview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



//! construct a row vector from given a diagview
template<typename eT>
inline Row<eT>::Row(const diagview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



//! construct a row vector from given a diagview
template<typename eT>
inline
const Row<eT>&
Row<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  
  //std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Row<eT>&
Row<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
Row<eT>::Row(const Op<T1, op_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Row<eT>&
Row<eT>::operator=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Row<eT>&
Row<eT>::operator*=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Row<eT>::Row(const Glue<T1, T2, glue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Row<eT>&
Row<eT>::operator=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Row<eT>&
Row<eT>::operator*=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
void
Row<eT>::set_size(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size(1,in_n_elem);
  }



template<typename eT>
inline
void
Row<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size( (std::min)( u32(1), in_n_rows), in_n_cols );
  
  arma_debug_check( (in_n_rows > 1), "Row::set_size(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros();
  }



template<typename eT>
inline
void
Row<eT>::zeros(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros(1,in_n_elem);
  }



template<typename eT>
inline
void
Row<eT>::zeros(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros( (std::min)( u32(1), in_n_rows), in_n_cols );
  arma_debug_check( (in_n_rows > 1), "Row<eT>::zeros(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Row<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::load(name,type);
  arma_debug_check( (Mat<eT>::n_rows > 1), "Row(): incompatible dimensions" );
  }



//! @}
