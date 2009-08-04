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


//! \addtogroup operator_minus
//! @{



//! unary -
template<typename T1>
arma_inline
const Op<T1, op_neg>
operator-
(const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1,op_neg>(X.get_ref());
  }



//! cancellation of two consecutive negations: -(-T1)
template<typename T1>
arma_inline
const T1&
operator-
(const Op<T1, op_neg>& X)
  {
  arma_extra_debug_sigprint();
  
  return X.m;
  }



//! Base - scalar
template<typename T1>
arma_inline
const Op<T1, op_scalar_minus_post>
operator-
(const Base<typename T1::elem_type,T1>& X, const typename T1::elem_type k)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_minus_post>(X.get_ref(), k);
  }



//! scalar - Base
template<typename T1>
arma_inline
const Op<T1, op_scalar_minus_pre>
operator-
(const typename T1::elem_type k, const Base<typename T1::elem_type,T1>& X)
  {
  arma_extra_debug_sigprint();
  
  return Op<T1, op_scalar_minus_pre>(X.get_ref(), k);
  }



//! Base - - Base = Base + Base
template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_plus>
operator-
  (
  const Base<typename T1::elem_type, T1                 >& X,
  const Base<typename T1::elem_type, Op<T2,op_neg> >& Y
  )
  {
  arma_extra_debug_sigprint();
  
  const T1& A = X.get_ref();
  const T2& B = (Y.get_ref()).m;
  
  return Glue<T1, T2, glue_plus>(A,B);
  }



//! Base - diagmat
template<typename T1, typename T2>
arma_inline
const Glue<T1, Op<T2,op_diagmat>, glue_minus_diag>
operator-
(const Base<typename T2::elem_type,T1>& X, const Op<T2,op_diagmat>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, Op<T1,op_diagmat>, glue_minus_diag>(X.get_ref(), Y);
  }



//! diagmat - Base
template<typename T1, typename T2>
arma_inline
const Glue< Op<T1,op_diagmat>, T2, glue_minus_diag>
operator-
(const Op<T1,op_diagmat>& X, const Base<typename T1::elem_type,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue< Op<T1,op_diagmat>, T2, glue_minus_diag>(X, Y.get_ref());
  }



//! Base - Op<T2,op_neg> = Base + T2
template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_plus>
operator-
(const Base<typename T2::elem_type,T1>& X, const Op<T2, op_neg>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_plus>(X.get_ref(), Y.m);
  }



//! diagmat - Op<T2,op_neg> = diagmat + T2
template<typename T1, typename T2>
arma_inline
const Glue< Op<T1,op_diagmat>, T2, glue_plus_diag>
operator-
  (
  const Base<typename T1::elem_type, Op<T1,op_diagmat> >& X,
  const Base<typename T1::elem_type, Op<T2,op_neg    > >& Y
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue< Op<T1,op_diagmat>, T2, glue_plus_diag>(X.get_ref(), (Y.get_ref()).m);
  }



//
// subtraction of Base objects with different element types
//



//! Base - Base
template<typename eT1, typename T1, typename eT2, typename T2>
arma_inline
Mat<typename promote_type<eT1,eT2>::result>
operator-
(const Base<eT1,T1>& X, const Base<eT2,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  promote_type<eT1,eT2>::check();
  
  const unwrap<T1> tmp1(X.get_ref());
  const unwrap<T2> tmp2(Y.get_ref());
  
  const Mat<eT1>& A = tmp1.M;
  const Mat<eT2>& B = tmp2.M;
  
  Mat< typename promote_type<eT1,eT2>::result > out;
  glue_minus::apply_mixed(out, A, B);
  
  return out;
  }



//
// subtraction of Base objects with same element types
//



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<std::complex<double>,T1>& X, const Base<std::complex<double>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<std::complex<float>,T1>& X, const Base<std::complex<float>,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<double,T1>& X, const Base<double,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<float,T1>& X, const Base<float,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<s32,T1>& X, const Base<s32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<u32,T1>& X, const Base<u32,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<s16,T1>& X, const Base<s16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<u16,T1>& X, const Base<u16,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<s8,T1>& X, const Base<s8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



template<typename T1, typename T2>
arma_inline
const Glue<T1, T2, glue_minus>
operator-
(const Base<u8,T1>& X, const Base<u8,T2>& Y)
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_minus>(X.get_ref(), Y.get_ref());
  }



//
// old operators
//


// //! unary -
// template<typename mT>
// arma_inline
// const Op<Mat<mT>, op_neg>
// operator-
// (const Mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Mat<mT>,op_neg>(X);
//   }
// 
// 
// 
// //! unary -
// template<typename T1, typename op_type>
// arma_inline
// const Op< Op<T1, op_type>, op_neg >
// operator-
// (const Op<T1, op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Op<T1, op_type>,op_neg>(X);
//   }
// 
// 
// 
// //! cancellation of two consecutive negations: -(-T1)
// template<typename T1>
// arma_inline
// const T1&
// operator-
// (const Op<T1, op_neg>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return X.m;
//   }
// 
// 
// 
// //! -(Glue<T1,T2,glue_type>)
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Op< Glue<T1,T2,glue_type>, op_neg >
// operator-
// (const Glue<T1,T2,glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op< Glue<T1,T2,glue_type>,op_neg >(X);
//   }
// 
// 
// 
// //! mat - scalar
// template<typename mT>
// arma_inline
// const Op<Mat<mT>, op_scalar_minus_post>
// operator-
// (const Mat<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Mat<mT>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! op - scalar
// template<typename T1, typename op_type>
// arma_inline
// const Op<Op<T1,op_type>, op_scalar_minus_post>
// operator-
// (const Op<T1,op_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Op<T1,op_type>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! op - scalar, level 2
// template<typename T1>
// arma_inline
// const Op<T1, op_scalar_minus_post>
// operator-
// (const Op<T1,op_scalar_plus>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<T1,op_scalar_minus_post>(X.m, X.aux - k);
//   }
// 
// 
// 
// //! glue - scalar
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Op<Glue<T1,T2,glue_type>, op_scalar_minus_post>
// operator-
// (const Glue<T1,T2,glue_type>& X, const typename T1::elem_type k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Glue<T1,T2,glue_type>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! scalar - mat
// template<typename mT>
// arma_inline
// const Op<Mat<mT>, op_scalar_minus_pre>
// operator-
// (const mT k, const Mat<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Mat<mT>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! scalar - op
// template<typename T1, typename op_type>
// arma_inline
// const Op<Op<T1,op_type>, op_scalar_minus_pre>
// operator-
// (const typename T1::elem_type k, const Op<T1,op_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Op<T1,op_type>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! scalar - glue
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Op<Glue<T1,T2,glue_type>, op_scalar_minus_pre>
// operator-
// (const typename T1::elem_type k, const Glue<T1,T2,glue_type>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<Glue<T1,T2,glue_type>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! mat - mat
// template<typename mT>
// arma_inline
// const Glue<Mat<mT>, Mat<mT>, glue_minus>
// operator-
// (const Mat<mT>& X, const Mat<mT>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<mT>, Mat<mT>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! mat - diagmat(T1)
// template<typename T1>
// arma_inline
// const Glue<Mat<typename T1::elem_type>, Op<T1,op_diagmat>, glue_minus_diag>
// operator-
// (const Mat<typename T1::elem_type>& X, const Op<T1,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T1::elem_type>, Op<T1,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - mat
// template<typename T1>
// arma_inline
// const Glue< Op<T1,op_diagmat>, Mat<typename T1::elem_type>, glue_minus_diag>
// operator-
// (const Op<T1,op_diagmat>& X, const Mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, Mat<typename T1::elem_type>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - diagmat(T2)
// template<typename T1, typename T2>
// arma_inline
// const Glue< Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_minus_diag>
// operator-
// (const Op<T1,op_diagmat>& X, const Op<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, Op<T2,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! mat - op
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Mat<typename T1::elem_type>, Op<T1,op_type>, glue_minus>
// operator-
// (const Mat<typename T1::elem_type>& X, const Op<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T1::elem_type>, Op<T1,op_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - op
// template<typename T1, typename T2, typename op_type>
// arma_inline
// const Glue< Op<T1,op_diagmat>, Op<T2,op_type>, glue_minus_diag>
// operator-
// (const Op<T1,op_diagmat>& X, const Op<T2,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, Op<T2,op_type>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! op - mat
// template<typename T1, typename op_type>
// arma_inline
// const Glue<Op<T1,op_type>, Mat<typename T1::elem_type>, glue_minus>
// operator-
// (const Op<T1,op_type>& X, const Mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Mat<typename T1::elem_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! op - diagmat(T2)
// template<typename T1, typename op_type, typename T2>
// arma_inline
// const Glue<Op<T1,op_type>, Op<T2,op_diagmat>, glue_minus_diag>
// operator-
// (const Op<T1,op_type>& X, const Op<T2,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Op<T2,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! op - glue
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// arma_inline
// const Glue<Op<T1,op_type>, Glue<T2, T3, glue_type>, glue_minus>
// operator-
// (const Op<T1,op_type>& X, const Glue<T2, T3, glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type>, Glue<T2, T3, glue_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! glue - op
// template<typename T1, typename op_type, typename T2, typename T3, typename glue_type>
// arma_inline
// const Glue<Glue<T2, T3, glue_type>, Op<T1,op_type>, glue_minus>
// operator-
// (const Glue<T2, T3, glue_type>& X, const Op<T1,op_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T2, T3, glue_type>, Op<T1,op_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! mat - glue
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue<Mat<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_minus>
// operator-
// (const Mat<typename T1::elem_type>& X, const Glue<T1, T2,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T1::elem_type>, Glue<T1, T2,glue_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! diagmat(T1) - glue
// template<typename T1, typename T2, typename T3, typename glue_type>
// arma_inline
// const Glue<Op<T1,op_diagmat>, Glue<T2,T3,glue_type>, glue_minus_diag>
// operator-
// (const Op<T1,op_diagmat>& X, const Glue<T2,T3,glue_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_diagmat>, Glue<T2,T3,glue_type>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! glue - mat
// template<typename T1, typename T2, typename glue_type>
// arma_inline
// const Glue< Glue<T1,T2,glue_type>, Mat<typename T1::elem_type>, glue_minus>
// operator-
// (const Glue<T1,T2,glue_type>& X, const Mat<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Glue<T1,T2,glue_type>, Mat<typename T1::elem_type>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! glue - diagmat(T1)
// template<typename T1, typename T2, typename glue_type, typename T3>
// arma_inline
// const Glue< Glue<T1,T2,glue_type>, Op<T3,op_diagmat>, glue_minus_diag>
// operator-
// (const Glue<T1,T2,glue_type>& X, const Op<T3,op_diagmat>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Glue<T1,T2,glue_type>, Op<T3,op_diagmat>, glue_minus_diag>(X,Y);
//   }
// 
// 
// 
// //! op - op
// template<typename T1, typename op_type1, typename T2, typename op_type2>
// arma_inline
// const Glue<Op<T1,op_type1>, Op<T2,op_type2>, glue_minus>
// operator-
// (const Op<T1,op_type1>& X, const Op<T2,op_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Op<T1,op_type1>, Op<T2,op_type2>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! glue - glue
// template<typename T1, typename T2, typename glue_type1, typename T3, typename T4, typename glue_type2>
// arma_inline
// const Glue<Glue<T1,T2,glue_type1>, Glue<T3,T4,glue_type2>, glue_minus>
// operator-
// (const Glue<T1,T2,glue_type1>& X, const Glue<T3,T4,glue_type2>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Glue<T1,T2,glue_type1>, Glue<T3,T4,glue_type2>, glue_minus>(X,Y);
//   }
// 
// 
// 
// //! mat - Op<T2,op_neg> = mat + T2
// template<typename T2>
// arma_inline
// const Glue<Mat<typename T2::elem_type>, T2, glue_plus>
// operator-
// (const Mat<typename T2::elem_type>& X, const Op<T2, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<Mat<typename T2::elem_type>, T2, glue_plus>(X,Y.m);
//   }
// 
// 
// 
// //! diagmat - Op<T2,op_neg> = diagmat + T2
// template<typename T1, typename T2>
// arma_inline
// const Glue< Op<T1,op_diagmat>, T2, glue_plus_diag>
// operator-
// (const Op<T1,op_diagmat>& X, const Op<T2, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1,op_diagmat>, T2, glue_plus_diag>(X,Y.m);
//   }
// 
// 
// 
// //! op - Op<T2,op_neg> = op + T2
// template<typename T1, typename op_type, typename T2>
// arma_inline
// const Glue<Op<T1, op_type>, T2, glue_plus>
// operator-
// (const Op<T1, op_type>& X, const Op<T2, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Op<T1, op_type>, T2, glue_plus>(X,Y.m);
//   }
// 
// 
// 
// //! glue - Op<T3,op_neg> = glue + T3
// template<typename T1, typename T2, typename glue_type, typename T3>
// arma_inline
// const Glue<Glue<T1,T2,glue_type>, T3, glue_plus>
// operator-
// (const Glue<T1,T2,glue_type>& X, const Op<T3, op_neg>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue< Glue<T1,T2,glue_type>, T3, glue_plus>(X,Y.m);
//   }
// 
// 
// 
// //! Base - subview
// template<typename T1>
// arma_inline
// const Glue<T1, subview<typename T1::elem_type>, glue_minus>
// operator-
// (const Base<T1>& X, const subview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<T1, subview<typename T1::elem_type>, glue_minus>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! subview - Base
// template<typename T1>
// arma_inline
// const Glue<subview<typename T1::elem_type>, T1, glue_minus>
// operator-
// (const subview<typename T1::elem_type>& X, const Base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<subview<typename T1::elem_type>, T1, glue_minus>(X,Y.get_ref());
//   }
// 
// 
// 
// //! Base - diagview
// template<typename T1>
// arma_inline
// const Glue<T1, diagview<typename T1::elem_type>, glue_minus>
// operator-
// (const Base<T1>& X, const diagview<typename T1::elem_type>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<T1, diagview<typename T1::elem_type>, glue_minus>(X.get_ref(),Y);
//   }
// 
// 
// 
// //! diagview - Base
// template<typename T1>
// arma_inline
// const Glue<diagview<typename T1::elem_type>, T1, glue_minus>
// operator-
// (const diagview<typename T1::elem_type>& X, const Base<T1>& Y)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Glue<diagview<typename T1::elem_type>, T1, glue_minus>(X,Y.get_ref());
//   }
// 
// 
// 
// //! scalar - subview
// template<typename mT>
// arma_inline
// const Op<subview<mT>, op_scalar_minus_pre>
// operator-
// (const mT k, const subview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<subview<mT>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! scalar - diagview
// template<typename mT>
// arma_inline
// const Op<diagview<mT>, op_scalar_minus_pre>
// operator-
// (const mT k, const diagview<mT>& X)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<diagview<mT>, op_scalar_minus_pre>(X,k);
//   }
// 
// 
// 
// //! subview - scalar
// template<typename mT>
// arma_inline
// const Op<subview<mT>, op_scalar_minus_post>
// operator-
// (const subview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<subview<mT>, op_scalar_minus_post>(X,k);
//   }
// 
// 
// 
// //! diagview - scalar
// template<typename mT>
// arma_inline
// const Op<diagview<mT>, op_scalar_minus_post>
// operator-
// (const diagview<mT>& X, const mT k)
//   {
//   arma_extra_debug_sigprint();
//   
//   return Op<diagview<mT>, op_scalar_minus_post>(X,k);
//   }



//! @}
