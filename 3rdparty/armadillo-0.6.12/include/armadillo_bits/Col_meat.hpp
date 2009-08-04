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


//! \addtogroup Col
//! @{


//! construct an empty column vector
template<typename eT>
inline
Col<eT>::Col()
  : Mat<eT>()
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector with the specified number of n_elem
template<typename eT>
inline
Col<eT>::Col(const u32 in_n_elem)
  : Mat<eT>(in_n_elem,1)
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector from specified text
template<typename eT>
inline
Col<eT>::Col(const char* text)
  : Mat<eT>(text)
  {
  arma_extra_debug_sigprint();
  
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from specified text
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(text);
  std::swap( access::rw(Mat<eT>::n_rows), access::rw(Mat<eT>::n_cols) );
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from a given column vector
template<typename eT>
inline
Col<eT>::Col(const Col<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  }



//! construct a column vector from a given column vector
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const Col<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  return *this;
  }



//! construct a column vector from a given matrix; the matrix must have exactly one column
template<typename eT>
inline
Col<eT>::Col(const Mat<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from a given matrix; the matrix must have exactly one column
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const Mat<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from a given auxillary array of eTs
template<typename eT>
inline
Col<eT>::Col(const eT* aux_mem, const u32 aux_length)
  {
  arma_extra_debug_sigprint();
  
  set_size(aux_length, 1);

  arma_check( (Mat<eT>::n_elem != aux_length), "Col::Col(): don't know how to handle the given array" );

  syslib::copy_elem( Mat<eT>::memptr(), aux_mem, Mat<eT>::n_elem );
  }



template<typename eT>
template<typename T1, typename T2>
inline
Col<eT>::Col
  (
  const Base<typename Col<eT>::pod_type, T1>& A,
  const Base<typename Col<eT>::pod_type, T2>& B
  )
  : Mat<eT>(A,B)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a submatrix; the submatrix must have exactly one column
template<typename eT>
inline
Col<eT>::Col(const subview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a submatrix; the submatrix must have exactly one column
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! construct a column vector from given a diagview
template<typename eT>
inline
Col<eT>::Col(const diagview<eT>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from given a diagview
template<typename eT>
inline
const Col<eT>&
Col<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  return *this;
  }



template<typename eT>
inline
const Col<eT>&
Col<eT>::operator*=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  return *this;
  }



//! construct a column vector from Op, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename op_type>
inline
Col<eT>::Col(const Op<T1, op_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from Op, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename op_type>
inline
const Col<eT>&
Col<eT>::operator=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col::operator=(): given matrix can't be interpreted as a column vector" );
  return *this;
  }



template<typename eT>
template<typename T1, typename op_type>
inline
const Col<eT>&
Col<eT>::operator*=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col::operator=(): incompatible dimensions" );
  return *this;
  }



//! construct a column vector from Glue, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Col<eT>::Col(const Glue<T1, T2, glue_type>& X)
  : Mat<eT>(X)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }



//! construct a column vector from Glue, i.e. run the previously delayed operations; the result of the operations must have exactly one column
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Col<eT>&
Col<eT>::operator=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Col<eT>&
Col<eT>::operator*=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::operator*=(X);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  
  return *this;
  }



//! change the number of n_rows
template<typename eT>
inline
void
Col<eT>::set_size(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size(in_n_elem,1);
  }



//! change the number of n_rows  (this function re-implements mat::set_size() in order to check the number of columns)
template<typename eT>
inline
void
Col<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::set_size( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  arma_debug_check( (in_n_cols > 1), "Col::set_size(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros();
  }



template<typename eT>
inline
void
Col<eT>::zeros(const u32 in_n_elem)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros(in_n_elem,1);
  }



template<typename eT>
inline
void
Col<eT>::zeros(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::zeros( in_n_rows, (std::min)( u32(1), in_n_cols ) );
  arma_debug_check( (in_n_cols > 1), "Col::zeros(): incompatible dimensions" );
  }



template<typename eT>
inline
void
Col<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>::load(name,type);
  arma_debug_check( (Mat<eT>::n_cols > 1), "Col(): incompatible dimensions" );
  }


//! @}
