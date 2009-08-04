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


//! \addtogroup operartor_times_dot
//! @{



//! Row vector multiplied by a column vector (i.e. a dot product)
template<typename eT>
arma_inline
eT
operator*
  (
  const Row<eT>& A,
  const Col<eT>& B
  )
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  
  }



//! Transpose of a column vector multiplied by a column vector (i.e. a dot product)
template<typename eT>
arma_inline
eT
operator*
  (
  const Op<Col<eT>, op_trans>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.m;
  const Col<eT>& B = Y;
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  
  }



//! Row vector multiplied by the transpose of a row vector (i.e. a dot product)
template<typename eT>
arma_inline
eT
operator*
  (
  const Row<eT>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Row<eT>& A = X;
  const Row<eT>& B = Y.m;
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  }
  

  
//! Transpose of a column vector multiplied by the transpose of a row vector (i.e. a dot product)
template<typename eT>
arma_inline
eT
operator*
  (
  const Op<Col<eT>, op_trans>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.m;
  const Row<eT>& B = Y.m;
  
  arma_debug_check( (A.n_elem != B.n_elem), "incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_elem, A.mem, B.mem);
  }



//
//
//



//! rowvec * mat * colvec 
template<typename eT>
inline
eT
operator*
  (
  const Glue<Row<eT>, Mat<eT>, glue_times_vec>& X,
  const Col<eT>& Y)
  {
  arma_extra_debug_sigprint();
  
  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * mat * colvec 
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>,op_trans>, Mat<eT>, glue_times>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! rowvec * mat * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue<Row<eT>, Mat<eT>, glue_times_vec>& X,
  const Op<Row<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * mat * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>,op_trans>, Mat<eT>, glue_times>& X,
  const Op<Row<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//
//
// 


//! rowvec * diagmat(mat) * colvec 
template<typename eT>
inline
eT
operator*
  (
  const Glue<Row<eT>, Op<Mat<eT>,op_diagmat>, glue_times_diag>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * diagmat(mat) * colvec 
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>,op_trans>, Op<Mat<eT>,op_diagmat>, glue_times_diag>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//! rowvec * diagmat(mat) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue<Row<eT>, Op<Mat<eT>,op_diagmat>, glue_times_diag>& X,
  const Op<Row<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * diagmat(mat) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>,op_trans>, Op<Mat<eT>,op_diagmat>, glue_times_diag>& X,
  const Op<Row<eT>,op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_diagmat_colvec(A.mem, B, C.mem);
  }



//
//
//

//! rowvec * diagmat(colvec or rowvec) * colvec
template<typename eT>
inline
eT
operator*
  (
  const Glue< Row<eT>, Op<Mat<eT>,op_diagmat_vec>, glue_times_vec>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//! trans(colvec) * diagmat(colvec or rowvec) * colvec
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>, op_trans>, Op<Mat<eT>,op_diagmat_vec>, glue_times>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//! rowvec * diagmat(colvec or rowvec) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Row<eT>, Op<Mat<eT>,op_diagmat_vec>, glue_times_vec>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//! trans(colvec) * diagmat(colvec or rowvec) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>, op_trans>, Op<Mat<eT>,op_diagmat_vec>, glue_times>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return op_dot::direct_dot(A.n_cols, A.mem, B.mem, C.mem);
  }



//
// 
// 



//! rowvec * inv(T1) * colvec
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const Glue< Row<typename T1::elem_type>, Op<T1, op_inv>, glue_times_vec>& X,
  const Col<typename T1::elem_type>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Row<eT>& A = X.A;
        Mat<eT>  B; op_inv::apply(B, X.B.m);
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(T1) * colvec
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const Glue< Op<Col<typename T1::elem_type>, op_trans>, Op<T1, op_inv>, glue_times>& X,
  const Col<typename T1::elem_type>& Y
  )
  {
  arma_extra_debug_sigprint();

  typedef typename T1::elem_type eT;
  
  const Col<eT>& A = X.A.m;
        Mat<eT>  B; op_inv::apply(B, X.B.m);
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! rowvec * inv(T1) * trans(rowvec)
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const Glue< Row<typename T1::elem_type>, Op<T1, op_inv>, glue_times_vec>& X,
  const Op<Row<typename T1::elem_type>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;

  const Row<eT>& A = X.A;
        Mat<eT>  B;  op_inv::apply(B, X.B.m);
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(T1) * trans(rowvec)
template<typename T1>
inline
typename T1::elem_type
operator*
  (
  const Glue< Op<Col<typename T1::elem_type>, op_trans>, Op<T1, op_inv>, glue_times>& X,
  const Op<Row<typename T1::elem_type>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Col<eT>& A = X.A.m;
        Mat<eT>  B;  op_inv::apply(B, X.B.m);
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_mat_colvec(A.mem, B, C.mem);
  }



// 
//
//



//! rowvec * inv(diagmat(mat)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const Glue< Row<eT>, Op< Op<Mat<eT>,op_diagmat>, op_inv>, glue_times_vec>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(mat)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>, op_trans>, Op< Op<Mat<eT>,op_diagmat>, op_inv>, glue_times>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//! rowvec * inv(diagmat(mat)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Row<eT>, Op< Op<Mat<eT>,op_diagmat>, op_inv>, glue_times_vec>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(mat)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>, op_trans>, Op< Op<Mat<eT>,op_diagmat>, op_inv>, glue_times>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_rows) || (B.n_rows != B.n_cols) || (B.n_cols != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagmat_colvec(A.mem, B, C.mem);
  }



//
//
//



//! rowvec * inv(diagmat(colvec or rowvec)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const Glue< Row<eT>, Op< Op<Mat<eT>,op_diagmat_vec>, op_inv>, glue_times_vec>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(colvec or rowvec)) * colvec
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>, op_trans>, Op< Op<Mat<eT>,op_diagmat_vec>, op_inv>, glue_times>& X,
  const Col<eT>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m.m;
  const Col<eT>& C = Y;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_rows), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//! rowvec * inv(diagmat(colvec or rowvec)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Row<eT>, Op< Op<Mat<eT>,op_diagmat_vec>, op_inv>, glue_times_vec>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();

  const Row<eT>& A = X.A;
  const Mat<eT>& B = X.B.m.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_cols != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//! trans(colvec) * inv(diagmat(colvec or rowvec)) * trans(rowvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue< Op<Col<eT>, op_trans>, Op< Op<Mat<eT>,op_diagmat_vec>, op_inv>, glue_times>& X,
  const Op<Row<eT>, op_trans>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A = X.A.m;
  const Mat<eT>& B = X.B.m.m;
  const Row<eT>& C = Y.m;
  
  arma_debug_check( (A.n_rows != B.n_elem) || (B.n_elem != C.n_cols), "operator*(): incompatible matrix dimensions" );
  
  return glue_times::direct_rowvec_invdiagvec_colvec(A.mem, B, C.mem);
  }



//
//
//



//! trans(colvec-colvec) * inv(diagmat(colvec or rowvec)) * (colvec-colvec)
template<typename eT>
inline
eT
operator*
  (
  const Glue
    <
    Op< Glue<Col<eT>, Col<eT>, glue_minus>, op_trans>,
    Op< Op<Mat<eT>,op_diagmat_vec>, op_inv>,
    glue_times
    >& X,
  const Glue<Col<eT>, Col<eT>, glue_minus>& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const Col<eT>& A1 = X.A.m.A;
  const Col<eT>& A2 = X.A.m.B;
  
  const Mat<eT>& B  = X.B.m.m;
  
  const Col<eT>& C1 = Y.A;
  const Col<eT>& C2 = Y.B;
  
  arma_debug_check
    (
    (A1.n_rows != A2.n_rows)
    ||
    (A2.n_rows != B.n_rows)
    ||
    (B.n_elem != C1.n_rows)
    ||
    (C1.n_rows != C2.n_rows)
    ,
    "operator*(): incompatible matrix dimensions"
    );
  
  
  if( (&A1 == &C1) && (&A2 == &C2) )
    {
    eT val = eT(0);
    for(u32 i=0; i<A1.n_rows; ++i)
      {
      const eT tmp = (A1.mem[i] - A2.mem[i]);
      val += (tmp*tmp) / B.mem[i];
      }
    return val;
    }
  else
    {
    eT val = eT(0);
    for(u32 i=0; i<A1.n_rows; ++i)
      {
      val += ( (A1.mem[i] - A2.mem[i]) * (C1.mem[i] - C2.mem[i]) ) / B.mem[i];
      }
    return val;
    }

  }



//! @}
