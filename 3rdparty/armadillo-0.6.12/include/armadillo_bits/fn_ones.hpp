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


//! \addtogroup fn_ones
//! @{



//! Delayed generation of a dense matrix with all elements set to one

inline
const Op<mat,op_ones_full>
ones(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Op<mat,op_ones_full>(n_rows, n_cols, 'j');
  }



template<typename mat_type>
inline
const Op<mat_type,op_ones_full>
ones(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Mat<mat_type>::value == false>::apply();
  
  return Op<mat_type,op_ones_full>(n_rows, n_cols, 'j');
  }



//! Generate a vector with all elements set to one
inline
const Op<colvec, op_ones_full>
ones(const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  return Op<colvec, op_ones_full>(n_elem, 1, 'j');
  }



template<typename vec_type>
inline
const Op<vec_type, op_ones_full>
ones(const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check< (is_Col<vec_type>::value == false) && (is_Row<vec_type>::value == false) >::apply();

  return Op<vec_type, op_ones_full>(n_elem, 0, 'j');
  }



//! Delayed generation of a diagonal matrix with the diagonal elements set to one
inline
const Op<mat,op_ones_diag>
eye(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  return Op<mat,op_ones_diag>(n_rows, n_cols, 'j');
  }



template<typename mat_type>
inline
const Op<mat_type,op_ones_diag>
eye(const u32 n_rows, const u32 n_cols)
  {
  arma_extra_debug_sigprint();
  
  arma_type_check<is_Mat_only<mat_type>::value == false>::apply();
  
  return Op<mat_type,op_ones_diag>(n_rows, n_cols, 'j');
  }


//! @}
