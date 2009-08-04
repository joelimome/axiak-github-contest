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




//! \addtogroup operator_times
//! @{

// "rowvec * colvec" and other combinations which result
// in a scalar are declared in "operator_times_dot.hpp"

//
// new operators

//! Base * scalar
template<typename T1>
arma_inline
const Op<T1, op_scalar_times>
operator*
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_times>(X.get_ref(),k);
  }



//! op * scalar, level 2
template<typename T1>
arma_inline
const Op<T1,op_scalar_times>
operator*
(const Op<T1,op_scalar_times>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_times>(X.m, X.aux * k);
  }



//! Op<mat,op_ones_full> * scalar
template<typename eT>
arma_inline
Mat<eT>
operator*
(const Op<Mat<eT>,op_ones_full>& X, const eT k)
  {
  arma_extra_debug_sigprint();
    
  Mat<eT> tmp(X.aux_u32_a, X.aux_u32_b);
  tmp.fill(k);
  
  return tmp;
  }



//! Op<mat,op_ones_diag> * scalar
template<typename eT>
arma_inline
Mat<eT>
operator*
(const Op<Mat<eT>,op_ones_diag>& X, const eT k)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.aux_u32_a, X.aux_u32_b);
  
  for(u32 i=0; i<out.n_rows; ++i)
    {
    out.at(i,i) = k;
    }
  
  return out;
  }



//! scalar * Base
template<typename T1>
arma_inline
const Op<T1, op_scalar_times>
operator*
(const typename T1::elem_type k, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_times>(X.get_ref(),k);  // NOTE: order is swapped
  }



//! scalar * Op<mat,op_ones_full>
template<typename eT>
arma_inline
Mat<eT>
operator*
(const eT k, const Op<Mat<eT>,op_ones_full>& X)
  {
  arma_extra_debug_sigprint();
    
  Mat<eT> tmp(X.aux_u32_a, X.aux_u32_b);
  tmp.fill(k);
  
  return tmp;
  }



//! scalar * Op<mat,op_ones_diag>
template<typename eT>
arma_inline
Mat<eT>
operator*
(const eT k, const Op<Mat<eT>,op_ones_diag>& X)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> out;
  out.zeros(X.aux_u32_a, X.aux_u32_b);
  
  for(u32 i=0; i<out.n_rows; ++i)
    out.at(i,i) = k;
  
  return out;
  }



//! Base * diagmat
template<typename T1, typename T2>
arma_inline
const Glue<T1, Op<T2,op_diagmat>, glue_times_diag>
operator*
(const Base<typename T2::elem_type, T1>& X, const Op<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T2,op_diagmat>, glue_times_diag>(X.get_ref(), Y);
  }



//! diagmat * Base
template<typename T1, typename T2>
arma_inline
const Glue<Op<T1,op_diagmat>, T2, glue_times_diag>
operator*
(const Op<T1,op_diagmat>& X, const Base<typename T1::elem_type, T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1,op_diagmat>, T2, glue_times_diag>(X, Y.get_ref());
  }



//! colvec * rowvec
template<typename eT>
arma_inline
const Glue<Col<eT>, Row<eT>, glue_times_vec>
operator*
(const Col<eT>& X, const Row<eT>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Col<eT>, Row<eT>, glue_times_vec>(X, Y);
  }



//! Base * colvec
template<typename T1>
arma_inline
const Glue<T1, Col<typename T1::elem_type>, glue_times_vec>
operator*
(const Base<typename T1::elem_type,T1>& X, const Col<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Col<typename T1::elem_type>, glue_times_vec>(X.get_ref(), Y);
  }



//! Base * rowvec
template<typename T1>
arma_inline
const Glue<T1, Row<typename T1::elem_type>, glue_times_vec>
operator*
(const Base<typename T1::elem_type,T1>& X, const Row<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Row<typename T1::elem_type>, glue_times_vec>(X.get_ref(), Y);
  }



//! colvec * Base
template<typename T1>
arma_inline
const Glue<Col<typename T1::elem_type>, T1, glue_times_vec>
operator*
(const Col<typename T1::elem_type>& X, const Base<typename T1::elem_type,T1>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Col<typename T1::elem_type>, T1, glue_times_vec>(X, Y.get_ref());
  }



//! rowvec * Base
template<typename T1>
arma_inline
const Glue<Row<typename T1::elem_type>, T1, glue_times_vec>
operator*
(const Row<typename T1::elem_type>& X, const Base<typename T1::elem_type,T1>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Row<typename T1::elem_type>, T1, glue_times_vec>(X, Y.get_ref());
  }



//! diagmat * colvec
template<typename T1>
arma_inline
const Glue<Op<T1, op_diagmat>, Col<typename T1::elem_type>, glue_times_diag>
operator*
(const Op<T1, op_diagmat>& X, const Col<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, Col<typename T1::elem_type>, glue_times_diag>(X, Y);
  }



//! diagmat * rowvec
template<typename T1>
arma_inline
const Glue<Op<T1, op_diagmat>, Row<typename T1::elem_type>, glue_times_diag>
operator*
(const Op<T1, op_diagmat>& X, const Row<typename T1::elem_type>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Op<T1, op_diagmat>, Row<typename T1::elem_type>, glue_times_diag>(X, Y);
  }



//! colvec * diagmat
template<typename T1>
arma_inline
const Glue<Col<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>
operator*
(const Col<typename T1::elem_type>& X, const Op<T1, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Col<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>(X, Y);
  }



//! rowvec * diagmat
template<typename T1>
arma_inline
const Glue<Row<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>
operator*
(const Row<typename T1::elem_type>& X, const Op<T1, op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<Row<typename T1::elem_type>, Op<T1, op_diagmat>, glue_times_diag>(X, Y);
  }



//
// multiplication of Base objects with different element types
//



//! Base * Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Mat<typename promote_type<eT1,eT2>::result>
operator*
(const Base<eT1,T1>& X, const Base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();

  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat< eT1 >& A = tmp1.M;
  const Mat< eT2 >& B = tmp2.M;
  
  Mat< typename promote_type<eT1,eT2>::result > out;
  
  glue_times::apply_mixed(out, A, B);
  
  return out;
  }



//
// multiplication of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<std::complex<double>,T1>& X, const Base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<std::complex<float>,T1>& X, const Base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<double,T1>& X, const Base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<float,T1>& X, const Base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<s32,T1>& X, const Base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<u32,T1>& X, const Base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<s16,T1>& X, const Base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<u16,T1>& X, const Base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<s8,T1>& X, const Base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_times>
operator*
(const Base<u8,T1>& X, const Base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_times>(X.get_ref(), Y.get_ref());
  }



// 
// old operators
//

// //! mat * scalar
// template<typename mT>
// arma_inline
// const Op<Mat<mT>, op_scalar_times>
// operator*
// (const Mat<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Mat<mT>, op_scalar_times>(X,k);
//   }
// 
// 
// 
// //! op * scalar
// template<typename T1, typename op_type>
// arma_inline
// const Op< Op<T1,op_type>, op_scalar_times>
// operator*
// (const Op<T1,op_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op< Op<T1,op_type>, op_scalar_times>(X,k);
//   }
// 
// 
// 
// //! op * scalar, level 2
// template<typename T1>
// arma_inline
// const Op<T1,op_scalar_times>
// operator*
// (const Op<T1,op_scalar_times>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<T1, op_scalar_times>(X.m, X.aux * k);
//   }
// 
// 
// 
// //! glue * scalar
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Op<Glue<T1,T2,glue_type>, op_scalar_times>
// operator*
// (const Glue<T1,T2,glue_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Glue<T1,T2,glue_type>, op_scalar_times>(X,k);
//   }
// 
// 
// 
// //! Op<mat,op_ones_full> * scalar
// template<typename mT>
// arma_inline
// Mat<mT>
// operator*
// (const Op<Mat<mT>,op_ones_full>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//     
//   Mat<mT> tmp(X.aux_u32_a, X.aux_u32_b);
//   tmp.fill(k);
//   
//   return tmp;
//   }
// 
// 
// 
// //! Op<mat,op_ones_diag> * scalar
// template<typename mT>
// arma_inline
// Mat<mT>
// operator*
// (const Op<Mat<mT>,op_ones_diag>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   Mat<mT> out;
//   out.zeros(X.aux_u32_a, X.aux_u32_b);
//   
//   for(u32 i=0; i<out.n_rows; ++i)
//     {
//     out.at(i,i) = k;
//     }
//   
//   return out;
//   }
// 
// 
// 
// //! scalar * mat
// template<typename mT>
// arma_inline
// const Op<Mat<mT>, op_scalar_times>
// operator*
// (const mT k, const Mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Mat<mT>, op_scalar_times>(X,k);  // NOTE: order is swapped
//   }
// 
// 
// 
// //! scalar * op
// template<typename T1, typename op_type>
// arma_inline
// const Op<Op<T1,op_type>, op_scalar_times>
// operator*
// (const typename T1::elem_type k, const Op<T1,op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Op<T1,op_type>, op_scalar_times>(X,k);  // NOTE: order is swapped
//   }
// 
// 
// 
// //! scalar * glue
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Op<Glue<T1,T2,glue_type>, op_scalar_times>
// operator*
// (const typename T1::elem_type k, const Glue<T1,T2,glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Glue<T1,T2,glue_type>, op_scalar_times>(X,k);  // NOTE: order is swapped
//   }
// 
// 
// 
// //! scalar * Op<mat,op_ones_full>
// template<typename mT>
// arma_inline
// Mat<mT>
// operator*
// (const mT k, const Op<Mat<mT>,op_ones_full>& X)
//   {
//   arma_extra_debug_sigprint();
//     
//   Mat<mT> tmp(X.aux_u32_a, X.aux_u32_b);
//   tmp.fill(k);
//   
//   return tmp;
//   }
// 
// 
// 
// //! scalar * Op<mat,op_ones_diag>
// template<typename mT>
// arma_inline
// Mat<mT>
// operator*
// (const mT k, const Op<Mat<mT>,op_ones_diag>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   Mat<mT> out;
//   out.zeros(X.aux_u32_a, X.aux_u32_b);
//   
//   for(u32 i=0; i<out.n_rows; ++i)
//     out.at(i,i) = k;
//   
//   return out;
//   }
// 
// 
// 
// //! mat * mat
// template<typename mT>
// arma_inline
// const Glue<Mat<mT>, Mat<mT>, glue_times>
// operator*
// (const Mat<mT>& X, const Mat<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<mT>, Mat<mT>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! mat * diagmat(T1)
// template<typename T1>
// arma_inline
// const Glue<Mat<typename T1::elem_type>, Op<T1,op_diagmat>, glue_times_diag>
// operator*
// (const Mat<typename T1::elem_type>& X, const Op<T1,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T1::elem_type>, Op<T1,op_diagmat>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) * mat
// template<typename T1>
// arma_inline
// const Glue<Op<T1,op_diagmat>, Mat<typename T1::elem_type>, glue_times_diag>
// operator*
// (const Op<T1,op_diagmat>& X, const Mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_diagmat>, Mat<typename T1::elem_type>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) * diagmat(T2)
// template<typename T1, typename T2>
// arma_inline
// const Glue< Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_times_diag>
// operator*
// (const Op<T1,op_diagmat>& X, const Op<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! mat * colvec
// template<typename mT>
// arma_inline
// const Glue<Mat<mT>, Col<mT>, glue_times>
// operator*
// (const Mat<mT>& X, const Col<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<mT>, Col<mT>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! mat * rowvec
// template<typename mT>
// arma_inline
// const Glue<Mat<mT>, Row<mT>, glue_times>
// operator*
// (const Mat<mT>& X, const Row<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<mT>, Row<mT>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! mat * op
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Mat<typename T1::elem_type>, Op<T1,op_type>, glue_times>
// operator*
// (const Mat<typename T1::elem_type>& X, const Op<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T1::elem_type>, Op<T1,op_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) * op
// template<typename T1, typename T2, typename op_type>
// arma_inline
// const Glue< Op<T1,op_diagmat>, Op<T2,op_type>, glue_times_diag>
// operator*
// (const Op<T1,op_diagmat>& X, const Op<T2,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, Op<T2,op_type>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! rowvec * mat
// template<typename mT>
// arma_inline
// const Glue<Row<mT>, Mat<mT>, glue_times>
// operator*
// (const Row<mT>& X, const Mat<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Row<mT>, Mat<mT>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! rowvec * op
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Row<typename T1::elem_type>, Op<T1,op_type>, glue_times>
// operator*
// (const Row<typename T1::elem_type>& X, const Op<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Row<typename T1::elem_type>, Op<T1,op_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! colvec * rowvec
// template<typename mT>
// arma_inline
// const Glue<Col<mT>, Row<mT>, glue_times>
// operator*
// (const Col<mT>& X, const Row<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Col<mT>, Row<mT>, glue_times>(X,Y);
//   }
// 
// 
// //! colvec * op
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Col<typename T1::elem_type>, Op<T1,op_type>, glue_times>
// operator*
// (const Col<typename T1::elem_type>& X, const Op<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Col<typename T1::elem_type>, Op<T1,op_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! op * mat
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Op<T1,op_type>, Mat<typename T1::elem_type>, glue_times>
// operator*
// (const Op<T1,op_type>& X, const Mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Mat<typename T1::elem_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! op * diagmat(T2)
// template<typename T1, typename op_type, typename T2>
// arma_inline
// const Glue<Op<T1,op_type>, Op<T2,op_diagmat>, glue_times_diag>
// operator*
// (const Op<T1,op_type>& X, const Op<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Op<T2,op_diagmat>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! op * colvec
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Op<T1,op_type>, Col<typename T1::elem_type>, glue_times>
// operator*
// (const Op<T1,op_type>& X, const Col<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Col<typename T1::elem_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! op * rowvec
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Op<T1,op_type>, Row<typename T1::elem_type>, glue_times>
// operator*
// (const Op<T1,op_type>& X, const Row<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Row<typename T1::elem_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! op * glue
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// arma_inline
// const Glue<Op<T1,op_type>, Glue<T2, T3, glue_type>, glue_times>
// operator*
// (const Op<T1,op_type>& X, const Glue<T2, T3, glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Glue<T2, T3, glue_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! glue * op
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// arma_inline
// const Glue<Glue<T2, T3, glue_type>, Op<T1,op_type>, glue_times>
// operator*
// (const Glue<T2, T3, glue_type>& X, const Op<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T2, T3, glue_type>, Op<T1,op_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! glue * mat
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Glue<T1, T2,glue_type>, Mat<typename T1::elem_type>, glue_times>
// operator*
// (const Glue<T1, T2,glue_type>& X, const Mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T1, T2,glue_type>, Mat<typename T1::elem_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! glue * diagmat(T3)
// template<typename T1, typename T2, typename glue_type, typename T3>
// arma_inline
// const Glue<Glue<T1, T2,glue_type>, Op<T3,op_diagmat>, glue_times_diag>
// operator*
// (const Glue<T1, T2,glue_type>& X, const Op<T3,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T1, T2,glue_type>, Op<T3,op_diagmat>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! glue * colvec
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Glue<T1, T2,glue_type>, Col<typename T1::elem_type>, glue_times>
// operator*
// (const Glue<T1, T2,glue_type>& X, const Col<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T1, T2,glue_type>, Col<typename T1::elem_type>, glue_times>(X,Y);
//   }
// 
// 
// //! glue * rowvec
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Glue<T1, T2,glue_type>, Row<typename T1::elem_type>, glue_times>
// operator*
// (const Glue<T1, T2,glue_type>& X, const Row<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T1, T2,glue_type>, Row<typename T1::elem_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! mat * glue
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Mat<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_times>
// operator*
// (const Mat<typename T1::elem_type>& X, const Glue<T1, T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) * glue
// template<typename T1, typename T2, typename T3, typename glue_type>
// arma_inline
// const Glue<Op<T1,op_diagmat>, Glue<T2,T3,glue_type>, glue_times_diag>
// operator*
// (const Op<T1,op_diagmat>& X, const Glue<T2,T3,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, Glue<T1, T2,glue_type>, glue_times_diag>(X,Y);
//   }
// 
// 
// 
// //! colvec * glue
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Col<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_times>
// operator*
// (const Col<typename T1::elem_type>& X, const Glue<T1, T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Col<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! rowvec * glue
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Row<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_times>
// operator*
// (const Row<typename T1::elem_type>& X, const Glue<T1, T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Row<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! op * op
// template<typename T1, typename op_type1, typename T2, typename op_type2>
// arma_inline
// const Glue<Op<T1,op_type1>, Op<T2,op_type2>, glue_times>
// operator*
// (const Op<T1,op_type1>& X, const Op<T2,op_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type1>, Op<T2,op_type2>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! glue * glue
// template<typename T1, typename T2, typename glue_type1, typename T3, typename T4, typename glue_type2>
// arma_inline
// const Glue<Glue<T1,T2,glue_type1>, Glue<T3,T4,glue_type2>, glue_times>
// operator*
// (const Glue<T1,T2,glue_type1>& X, const Glue<T3,T4,glue_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T1,T2,glue_type1>, Glue<T3,T4,glue_type2>, glue_times>(X,Y);
//   }
// 
// 
// 
// //! Base * subview
// template<typename T1>
// arma_inline
// const Glue<T1, subview<typename T1::elem_type>, glue_times>
// operator*
// (const Base<T1>& X, const subview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<T1, subview<typename T1::elem_type>, glue_times>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! subview * Base
// template<typename T1>
// arma_inline
// const Glue<subview<typename T1::elem_type>, T1, glue_times>
// operator*
// (const subview<typename T1::elem_type>& X, const Base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<subview<typename T1::elem_type>, T1, glue_times>(X,Y.get_ref());
//   }
// 
// 
// 
// //! Base * diagview
// template<typename T1>
// arma_inline
// const Glue<T1, diagview<typename T1::elem_type>, glue_times>
// operator*
// (const Base<T1>& X, const diagview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<T1, diagview<typename T1::elem_type>, glue_times>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! diagview * Base
// template<typename T1>
// arma_inline
// const Glue<diagview<typename T1::elem_type>, T1, glue_times>
// operator*
// (const diagview<typename T1::elem_type>& X, const Base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<diagview<typename T1::elem_type>, T1, glue_times>(X,Y.get_ref());
//   }
// 
// 
// 
// //! scalar * subview
// template<typename mT>
// arma_inline
// const Op<subview<mT>, op_scalar_times>
// operator*
// (const mT k, const subview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<subview<mT>, op_scalar_times>(X,k);
//   }
// 
// 
// 
// //! scalar * diagview
// template<typename mT>
// arma_inline
// const Op<diagview<mT>, op_scalar_times>
// operator*
// (const mT k, const diagview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<diagview<mT>, op_scalar_times>(X,k);
//   }
// 
// 
// 
// //! subview * scalar
// template<typename mT>
// arma_inline
// const Op<subview<mT>, op_scalar_times>
// operator*
// (const subview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<subview<mT>, op_scalar_times>(X,k);
//   }
// 
// 
// 
// //! diagview * scalar
// template<typename mT>
// arma_inline
// const Op<diagview<mT>, op_scalar_times>
// operator*
// (const diagview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<diagview<mT>, op_scalar_times>(X,k);
//   }
// 


//! @}
