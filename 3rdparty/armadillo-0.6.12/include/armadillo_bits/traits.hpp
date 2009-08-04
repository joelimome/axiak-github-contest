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


//! \addtogroup traits
//! @{


template<typename T1>
struct get_pod_type
  { typedef T1 pod_type; };

template<typename T2>
struct get_pod_type< std::complex<T2> >
  { typedef T2 pod_type; };



template<typename T>
struct is_Mat_only
  { static const bool value = false; };

template<typename eT>
struct is_Mat_only< Mat<eT> >
  { static const bool value = true; };



template<typename T>
struct is_Mat
  { static const bool value = false; };

template<typename eT>
struct is_Mat< Mat<eT> >
  { static const bool value = true; };

template<typename eT>
struct is_Mat< Row<eT> >
  { static const bool value = true; };

template<typename eT>
struct is_Mat< Col<eT> >
  { static const bool value = true; };



template<typename T>
struct is_Row
  { static const bool value = false; };

template<typename eT>
struct is_Row< Row<eT> >
  { static const bool value = true; };



template<typename T>
struct is_Col
  { static const bool value = false; };

template<typename eT>
struct is_Col< Col<eT> >
  { static const bool value = true; };






template<typename T>
struct is_subview
  { static const bool value = false; };

template<typename eT>
struct is_subview< subview<eT> >
  { static const bool value = true; };


template<typename T>
struct is_diagview
  { static const bool value = false; };

template<typename eT>
struct is_diagview< diagview<eT> >
  { static const bool value = true; };


//
//
//


template<typename T>
struct is_Op
  { static const bool value = false; };
 
template<typename T1, typename op_type>
struct is_Op< Op<T1,op_type> >
  { static const bool value = true; };
 
 
template<typename T>
struct is_Glue
  { static const bool value = false; };
 
template<typename T1, typename T2, typename glue_type>
struct is_Glue< Glue<T1,T2,glue_type> >
  { static const bool value = true; };


template<typename T>
struct is_glue_times
  { static const bool value = false; };

template<typename T1, typename T2>
struct is_glue_times< Glue<T1,T2,glue_times> >
  { static const bool value = true; };


//
//
//


template<typename T1>
struct is_arma_type
  {
  static const bool value
  =  is_Mat<T1>::value
  || is_Op<T1>::value
  || is_Glue<T1>::value
  || is_subview<T1>::value
  || is_diagview<T1>::value
  ;
  };



//
//
//


template<typename T1, typename T2>
struct is_same_type
  { static const bool value = false; };


template<typename T1>
struct is_same_type<T1,T1>
  { static const bool value = true; };



template<typename T1, typename T2>
struct isnt_same_type
  {
  static const bool value = true;
  
  inline static void check()
    {
    arma_static_assert<false> ERROR___TYPE_MISMATCH;
    ERROR___TYPE_MISMATCH = ERROR___TYPE_MISMATCH;
    }
  };


template<typename T1>
struct isnt_same_type<T1,T1>
  {
  static const bool value = false;
  
  arma_inline static void check() {}
  };


//
//
//


template<typename T1>
struct isnt_fltpt
  {
  static const bool value = true;
  
  inline static void check()
    {
    arma_static_assert<false> ERROR___TYPE_MISMATCH;
    ERROR___TYPE_MISMATCH = ERROR___TYPE_MISMATCH;
    }
  };



struct isnt_fltpt_false
  {
  static const bool value = false;
  
  arma_inline static void check() {}
  };



template<> struct isnt_fltpt< float >                : public isnt_fltpt_false {};
template<> struct isnt_fltpt< double >               : public isnt_fltpt_false {};
template<> struct isnt_fltpt< long double >          : public isnt_fltpt_false {};
template<> struct isnt_fltpt< std::complex<float> >  : public isnt_fltpt_false {};
template<> struct isnt_fltpt< std::complex<double> > : public isnt_fltpt_false {};



template<typename T1>
struct is_u8
  { static const bool value = false; };

template<>
struct is_u8<u8>
  { static const bool value = true; };



template<typename T1>
struct is_s8
  { static const bool value = false; };

template<>
struct is_s8<s8>
  { static const bool value = true; };



template<typename T1>
struct is_u16
  { static const bool value = false; };

template<>
struct is_u16<u16>
  { static const bool value = true; };



template<typename T1>
struct is_s16
  { static const bool value = false; };

template<>
struct is_s16<s16>
  { static const bool value = true; };



template<typename T1>
struct is_u32
  { static const bool value = false; };

template<>
struct is_u32<u32>
  { static const bool value = true; };



template<typename T1>
struct is_s32
  { static const bool value = false; };

template<>
struct is_s32<s32>
  { static const bool value = true; };




template<typename T1>
struct is_float
  { static const bool value = false; };

template<>
struct is_float<float>
  { static const bool value = true; };



template<typename T1>
struct is_double
  { static const bool value = false; };

template<>
struct is_double<double>
  { static const bool value = true; };



template<typename T1>
struct is_complex
  { static const bool value = false; };

// template<>
template<typename eT>
struct is_complex< std::complex<eT> >
  { static const bool value = true; };



template<typename T1>
struct is_complex_float
  { static const bool value = false; };

template<>
struct is_complex< std::complex<float> >
  { static const bool value = true; };



template<typename T1>
struct is_complex_double
  { static const bool value = false; };

template<>
struct is_complex_double< std::complex<double> >
  { static const bool value = true; };




//! check for a weird implementation of the std::complex class
template<typename T1>
struct is_supported_complex
  { static const bool value = false; };

//template<>
template<typename eT>
struct is_supported_complex< std::complex<eT> >
  { static const bool value = ( sizeof(std::complex<eT>) == 2*sizeof(eT) ); };



template<typename T1>
struct is_supported_complex_float
  { static const bool value = false; };

template<>
struct is_supported_complex_float< std::complex<float> >
  { static const bool value = ( sizeof(std::complex<float>) == 2*sizeof(float) ); };



template<typename T1>
struct is_supported_complex_double
  { static const bool value = false; };

template<>
struct is_supported_complex_double< std::complex<double> >
  { static const bool value = ( sizeof(std::complex<double>) == 2*sizeof(double) ); };



template<typename T1>
struct is_supported_elem_type
  {
  static const bool value = \
    is_u8<T1>::value ||
    is_s8<T1>::value ||
    is_u16<T1>::value ||
    is_s16<T1>::value ||
    is_u32<T1>::value ||
    is_s32<T1>::value ||
    is_float<T1>::value ||
    is_double<T1>::value ||
    is_supported_complex_float<T1>::value ||
    is_supported_complex_double<T1>::value;
  };



template<typename T1>
struct isnt_supported_elem_type
  {
  static const bool value = true;
  
  inline static void check()
    {
    arma_static_assert<false> ERROR___UNSUPPORTED_ELEMENT_TYPE;
    ERROR___UNSUPPORTED_ELEMENT_TYPE = ERROR___UNSUPPORTED_ELEMENT_TYPE;
    }
  };



struct isnt_supported_elem_type_false
  {
  static const bool value = false;
  
  arma_inline static void check() {}
  };



template<> struct isnt_supported_elem_type< u8 >                   : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< s8 >                   : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< u16 >                  : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< s16 >                  : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< u32 >                  : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< s32 >                  : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< float >                : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< double >               : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< std::complex<float> >  : public isnt_supported_elem_type_false {};
template<> struct isnt_supported_elem_type< std::complex<double> > : public isnt_supported_elem_type_false {};



template<typename T1>
struct is_supported_blas_type
  {
  static const bool value = \
    is_float<T1>::value ||
    is_double<T1>::value ||
    is_supported_complex_float<T1>::value ||
    is_supported_complex_double<T1>::value;
  };



template<typename T>
struct is_signed
  {
  static const bool value = (T(-1) < T(0));
  };



template<typename T>
struct is_non_integral
  {
  static const bool value = (T(1.0) != T(1.1));
  };



//
// type promotion

template<typename T1, typename T2>
struct promote_type
  {
  inline static void check()
    {
    arma_static_assert<false> ERROR___UNSUPPORTED_MIXTURE_OF_TYPES;
    ERROR___UNSUPPORTED_MIXTURE_OF_TYPES = ERROR___UNSUPPORTED_MIXTURE_OF_TYPES;
    }
  
  typedef T1 result;
  };



struct promote_type_ok
  {
  arma_inline static void check() {}
  };


template<typename T>
struct promote_type<T, T> : public promote_type_ok { typedef T result; };

template<typename T>
struct promote_type<std::complex<T>,     T>      : public promote_type_ok { typedef std::complex<T> result; };

template<>
struct promote_type<std::complex<double>, std::complex<float> > : public promote_type_ok { typedef std::complex<double> result; };

template<>
struct promote_type<std::complex<double>, float> : public promote_type_ok { typedef std::complex<double> result; };

template<>
struct promote_type<std::complex<float>, double> : public promote_type_ok { typedef std::complex<double> result; };

template<typename T>
struct promote_type<std::complex<T>, s32> : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<std::complex<T>, u32> : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<std::complex<T>, s16> : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<std::complex<T>, u16> : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<std::complex<T>, s8>  : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<std::complex<T>, u8>  : public promote_type_ok { typedef std::complex<T> result; };


template<> struct promote_type<double, float> : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, s32  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, u32  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, s16  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, u16  > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, s8   > : public promote_type_ok { typedef double result; };
template<> struct promote_type<double, u8   > : public promote_type_ok { typedef double result; };

template<> struct promote_type<float, s32> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, u32> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, s16> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, u16> : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, s8 > : public promote_type_ok { typedef float result; };
template<> struct promote_type<float, u8 > : public promote_type_ok { typedef float result; };

template<> struct promote_type<s32, u32> : public promote_type_ok { typedef s32 result; };  // float ?  
template<> struct promote_type<s32, s16> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s32, u16> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s32, s8 > : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s32, u8 > : public promote_type_ok { typedef s32 result; };

template<> struct promote_type<u32, s16> : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u32, u16> : public promote_type_ok { typedef u32 result; };
template<> struct promote_type<u32, s8 > : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u32, u8 > : public promote_type_ok { typedef u32 result; };

template<> struct promote_type<s16, u16> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<s16, s8 > : public promote_type_ok { typedef s16 result; };
template<> struct promote_type<s16, u8 > : public promote_type_ok { typedef s16 result; };

template<> struct promote_type<u16, s8> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<u16, u8> : public promote_type_ok { typedef u16 result; };

template<> struct promote_type<s8, u8> : public promote_type_ok { typedef s8 result; };  // s16 ?




//
// type promotion, mirrored versions

template<typename T>
struct promote_type<T, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };

template<>
struct promote_type<std::complex<float>, std::complex<double> > : public promote_type_ok { typedef std::complex<double> result; };

template<>
struct promote_type<float, std::complex<double> > : public promote_type_ok { typedef std::complex<double> result; };

template<>
struct promote_type<double, std::complex<float> > : public promote_type_ok { typedef std::complex<double> result; };

template<typename T>
struct promote_type<s32, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<u32, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<s16, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<u16, std::complex<T> > : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<s8, std::complex<T> >  : public promote_type_ok { typedef std::complex<T> result; };

template<typename T>
struct promote_type<u8, std::complex<T> >  : public promote_type_ok { typedef std::complex<T> result; };


template<> struct promote_type<float, double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<s32  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<u32  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<s16  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<u16  , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<s8   , double> : public promote_type_ok { typedef double result; };
template<> struct promote_type<u8   , double> : public promote_type_ok { typedef double result; };

template<> struct promote_type<s32, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<u32, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<s16, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<u16, float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<s8 , float> : public promote_type_ok { typedef float result; };
template<> struct promote_type<u8 , float> : public promote_type_ok { typedef float result; };

template<> struct promote_type<u32, s32> : public promote_type_ok { typedef s32 result; };  // float ?  
template<> struct promote_type<s16, s32> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<u16, s32> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<s8 , s32> : public promote_type_ok { typedef s32 result; };
template<> struct promote_type<u8 , s32> : public promote_type_ok { typedef s32 result; };

template<> struct promote_type<s16, u32> : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u16, u32> : public promote_type_ok { typedef u32 result; };
template<> struct promote_type<s8 , u32> : public promote_type_ok { typedef s32 result; };  // float ?
template<> struct promote_type<u8 , u32> : public promote_type_ok { typedef u32 result; };

template<> struct promote_type<u16, s16> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<s8 , s16> : public promote_type_ok { typedef s16 result; };
template<> struct promote_type<u8 , s16> : public promote_type_ok { typedef s16 result; };

template<> struct promote_type<s8, u16> : public promote_type_ok { typedef s16 result; };  // s32 ?
template<> struct promote_type<u8, u16> : public promote_type_ok { typedef u16 result; };

template<> struct promote_type<u8, s8> : public promote_type_ok { typedef s8 result; };  // s16 ?


//
//
//


//! upgrade_val is used to ensure an operation such as multiplication is possible between two types.
//! values are upgraded only where necessary.
//! used by:
//! glue_div::apply_mixed(), glue_minus::apply_mixed(), glue_plus::apply_mixed(), glue_schur::apply_mixed()
//! and glue_times::apply_mixed() via gemm_mixed().

template<typename T1, typename T2>
struct upgrade_val
  {
  typedef typename promote_type<T1,T2>::result T1_result;
  typedef typename promote_type<T1,T2>::result T2_result;
  
  arma_inline
  static
  const typename promote_type<T1,T2>::result
  apply(const T1 x)
    {
    typedef typename promote_type<T1,T2>::result out_type;
    return out_type(x);
    }
  
  arma_inline
  static
  const typename promote_type<T1,T2>::result
  apply(const T2 x)
    {
    typedef typename promote_type<T1,T2>::result out_type;
    return out_type(x);
    }
  
  };


// template<>
template<typename T>
struct upgrade_val<T,T>
  {
  typedef T T1_result;
  typedef T T2_result;
  
  arma_inline static const T& apply(const T& x) { return x; }
  };


//! upgrade a type to allow multiplication with a complex type
//! e.g. the int in "int * complex<double>" is upgraded to a double
// template<>
template<typename T, typename T2>
struct upgrade_val< std::complex<T>, T2 >
  {
  typedef std::complex<T> T1_result;
  typedef T               T2_result;
  
  arma_inline static const std::complex<T>& apply(const std::complex<T>& x) { return x;    }
  arma_inline static const T                apply(const T2 x)               { return T(x); }
  };


// template<>
template<typename T1, typename T>
struct upgrade_val< T1, std::complex<T> >
  {
  typedef T               T1_result;
  typedef std::complex<T> T2_result;
  
  arma_inline static const T                apply(const T1 x)               { return T(x); }
  arma_inline static const std::complex<T>& apply(const std::complex<T>& x) { return x;    }
  };


//! ensure we don't lose precision when multiplying a complex number with a higher precision real number
template<>
struct upgrade_val< std::complex<float>, double >
  {
  typedef std::complex<double> T1_result;
  typedef double               T2_result;
  
  arma_inline static const std::complex<double> apply(const std::complex<float>& x) { return std::complex<double>(x); }
  arma_inline static const double               apply(const double x)               { return x; }
  };


template<>
struct upgrade_val< double, std::complex<float> >
  {
  typedef double              T1_result;
  typedef std::complex<float> T2_result;
  
  arma_inline static const double               apply(const double x)               { return x; }
  arma_inline static const std::complex<double> apply(const std::complex<float>& x) { return std::complex<double>(x); }
  };


//! ensure we don't lose precision when multiplying complex numbers with different underlying types
template<>
struct upgrade_val< std::complex<float>, std::complex<double> >
  {
  typedef std::complex<double> T1_result;
  typedef std::complex<double> T2_result;
  
  arma_inline static const std::complex<double>  apply(const std::complex<float>&  x) { return std::complex<double>(x); }
  arma_inline static const std::complex<double>& apply(const std::complex<double>& x) { return x; }
  };


template<>
struct upgrade_val< std::complex<double>, std::complex<float> >
  {
  typedef std::complex<double> T1_result;
  typedef std::complex<double> T2_result;
  
  arma_inline static const std::complex<double>& apply(const std::complex<double>& x) { return x; }
  arma_inline static const std::complex<double>  apply(const std::complex<float>&  x) { return std::complex<double>(x); }
  };


//! work around limitations in the complex class (at least as present in gcc 4.1 & 4.3)
template<>
struct upgrade_val< std::complex<double>, float >
  {
  typedef std::complex<double> T1_result;
  typedef double               T2_result;
  
  arma_inline static const std::complex<double>& apply(const std::complex<double>& x) { return x; }
  arma_inline static const double                apply(const float x)                 { return double(x); }
  };


template<>
struct upgrade_val< float, std::complex<double> >
  {
  typedef double               T1_result;
  typedef std::complex<double> T2_result;
  
  arma_inline static const double                apply(const float x)                 { return double(x); }
  arma_inline static const std::complex<double>& apply(const std::complex<double>& x) { return x; }
  };


//! @}
