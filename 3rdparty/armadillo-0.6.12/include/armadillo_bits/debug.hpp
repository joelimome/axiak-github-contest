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


//! \addtogroup debug
//! @{


//
// arma_print


inline
void
arma_print()
  {
  std::cout << std::endl;
  }


template<typename T1>
inline
void
arma_print(const T1& x)
  {
  std::cout << x << std::endl;
  }



template<typename T1, typename T2>
inline
void
arma_print(const T1& x, const T2& y)
  {
  std::cout << x << y << std::endl;
  }



#ifdef ARMA_USE_BOOST
  template<typename T1>
  inline
  void
  arma_print(const arma_boost::basic_format<T1>& x)
    {
    std::cout << x << std::endl;
    }
#else
  template<typename T1, typename T2>
  inline
  void
  arma_print(const arma_boost::basic_format<T1,T2>& x)
    {
    std::cout << x << std::endl;
    }
#endif



//
// arma_sigprint

//! print a message on cout, with a preceding @ character.
//! used for printing the signature of a function
//! (see the arma_extra_debug_sigprint macro) 
inline
void
arma_sigprint(const char* x)
  {
  std::cout << "@ " << x;
  }



//
// arma_bktprint


inline
void
arma_bktprint()
  {
  std::cout << std::endl;
  }


template<typename T1>
inline
void
arma_bktprint(const T1& x)
  {
  std::cout << " [" << x << "]" << std::endl;
  }



template<typename T1, typename T2>
inline
void
arma_bktprint(const T1& x, const T2& y)
  {
  std::cout << " [" << x << y << "]" << std::endl;
  }



#ifdef ARMA_USE_BOOST
  template<typename T1>
  inline
  void
  arma_bktprint(const arma_boost::basic_format<T1>& x)
    {
    std::cout << " [" << x << "]" << std::endl;
    }
#else
  template<typename T1, typename T2>
  inline
  void
  arma_bktprint(const arma_boost::basic_format<T1,T2>& x)
    {
    std::cout << " [" << x << "]" << std::endl;
    }
#endif



//
// arma_thisprint


inline
void
arma_thisprint(void* this_ptr)
  {
  std::cout << " [this = " << this_ptr << "]" << std::endl;
  }



//
// arma_warn

//! if state is true, print a message on cout
template<typename T1>
inline
void
arma_hot
arma_warn(const bool state, const T1& x)
  {
  if(state==true)
    {
    arma_print(x);
    }
  }


template<typename T1, typename T2>
inline
void
arma_hot
arma_warn(const bool state, const T1& x, const T2& y)
  {
  if(state==true)
    {
    arma_print(x,y);
    }
  }


#ifdef ARMA_USE_BOOST
  template<typename T1>
  inline
  void
  arma_hot
  arma_warn(const bool state, const arma_boost::basic_format<T1>& x)
    {
    if(state==true)
      arma_print(x);
    }
#else
  template<typename T1, typename T2>
  inline
  void
  arma_hot
  arma_warn(const bool state, const arma_boost::basic_format<T1,T2>& x)
    {
    if(state==true)
      arma_print(x);
    }
#endif



//
// arma_check

//! if state is true, throw a run-time error exception
template<typename T1>
inline
void
arma_hot
arma_check(const bool state, const T1& x)
  {
  if(state==true)
    {
    throw std::runtime_error(x);
    }
  }


template<typename T1, typename T2>
inline
void
arma_hot
arma_check(const bool state, const T1& x, const T2& y)
  {
  if(state==true)
    {
    throw std::runtime_error( std::string(x) + std::string(y) );
    }
  }


#ifdef ARMA_USE_BOOST
  template<typename T1>
  inline
  void
  arma_hot
  arma_check(const bool state, const arma_boost::basic_format<T1>& x)
    {
    if(state==true)
      {
      throw std::runtime_error(str(x));
      }
    }
#else
  template<typename T1, typename T2>
  inline
  void
  arma_hot
  arma_check(const bool state, const arma_boost::basic_format<T1,T2>& x)
    {
    if(state==true)
      {
      throw std::runtime_error(str(x));
      }
    }
#endif



inline
std::string
arma_incompat_size_string(const u32 A_n_rows, const u32 A_n_cols, const u32 B_n_rows, const u32 B_n_cols, const char* x)
  {
  std::stringstream tmp;
  
  tmp << x << ": incompatible matrix dimensions: (" << A_n_rows << "," << A_n_cols << ") and (" << B_n_rows << "," << B_n_cols << ")";
  
  return tmp.str();
  }



inline
void
arma_hot
arma_assert_same_size(const u32 A_n_rows, const u32 A_n_cols, const u32 B_n_rows, const u32 B_n_cols, const char* x)
  {
  if( (A_n_rows != B_n_rows) || (A_n_cols != B_n_cols) )
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x)
      );
    }
  }



//! if given matrices have different sizes, throw a run-time error exception
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const Mat<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_same_size(const subview<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if( (A.n_rows != B.n_rows) || (A.n_cols != B.n_cols) )
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



inline
void
arma_hot
arma_assert_mul_size(const u32 A_n_rows, const u32 A_n_cols, const u32 B_n_rows, const u32 B_n_cols, const char* x)
  {
  if(A_n_cols != B_n_rows)
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A_n_rows, A_n_cols, B_n_rows, B_n_cols, x)
      );
    }
  }



//! if given matrices are incompatible for multiplication, throw a run-time error exception
template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const Mat<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const Mat<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const subview<eT1>& A, const Mat<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



template<typename eT1, typename eT2>
inline
void
arma_hot
arma_assert_mul_size(const subview<eT1>& A, const subview<eT2>& B, const char* x)
  {
  if(A.n_cols != B.n_rows)
    {
    throw std::runtime_error
      (
      arma_incompat_size_string(A.n_rows, A.n_cols, B.n_rows, B.n_cols, x)
      );
    }
  }



//
// arma_stop

//! throw a run-time error exception
template<typename T1>
inline
void
arma_stop(const T1& x)
  {
  arma_check(true, x);
  }



//
// macros


#define ARMA_STRING1(x) #x
#define ARMA_STRING2(x) ARMA_STRING1(x)
#define ARMA_FILELINE  __FILE__ ": " ARMA_STRING2(__LINE__)


#if defined (__GNUG__)
  #define ARMA_FNSIG  __PRETTY_FUNCTION__
#elif defined (_MSC_VER)
  #define ARMA_FNSIG  __FUNCSIG__ 
#elif defined (ARMA_USE_BOOST)
  #define ARMA_FNSIG  BOOST_CURRENT_FUNCTION  
#else 
  #define ARMA_FNSIG  "(unknown)"
#endif



#if !defined(ARMA_NO_DEBUG) && !defined(NDEBUG)
  
  #define arma_debug_print            arma_print
  #define arma_debug_warn             arma_warn
  #define arma_debug_check            arma_check
  #define arma_debug_assert_same_size arma_assert_same_size
  #define arma_debug_assert_mul_size  arma_assert_mul_size
  
#else
  
  #undef ARMA_EXTRA_DEBUG
  
  #define arma_debug_print            true ? (void)0 : arma_print
  #define arma_debug_warn             true ? (void)0 : arma_warn
  #define arma_debug_check            true ? (void)0 : arma_check
  #define arma_debug_assert_same_size true ? (void)0 : arma_assert_same_size
  #define arma_debug_assert_mul_size  true ? (void)0 : arma_assert_mul_size

#endif


#if defined(ARMA_EXTRA_DEBUG)
  
  #define arma_extra_debug_sigprint       arma_sigprint(ARMA_FNSIG); arma_bktprint
  #define arma_extra_debug_sigprint_this  arma_sigprint(ARMA_FNSIG); arma_thisprint
  #define arma_extra_debug_print          arma_print
  #define arma_extra_debug_warn           arma_warn
  #define arma_extra_debug_check          arma_check

#else
  
  #define arma_extra_debug_sigprint        true ? (void)0 : arma_bktprint
  #define arma_extra_debug_sigprint_this   true ? (void)0 : arma_thisprint
  #define arma_extra_debug_print           true ? (void)0 : arma_print
  #define arma_extra_debug_warn            true ? (void)0 : arma_warn
  #define arma_extra_debug_check           true ? (void)0 : arma_check
 
#endif




#if defined(ARMA_EXTRA_DEBUG)

  namespace junk
    {
    class arma_first_extra_debug_message
      {
      public:
      
      inline
      arma_first_extra_debug_message()
        {
        std::cout << "@ ---" << '\n';
        std::cout << "@ Armadillo " << arma_version::major << '.' << arma_version::minor << '.' << arma_version::patch << '\n';
        std::cout << "@ arma_config::atlas      = " << arma_config::atlas      << '\n';
        std::cout << "@ arma_config::lapack     = " << arma_config::lapack     << '\n';
        std::cout << "@ arma_config::blas       = " << arma_config::blas       << '\n';
        std::cout << "@ arma_config::boost      = " << arma_config::boost      << '\n';
        std::cout << "@ arma_config::boost_date = " << arma_config::boost_date << '\n';
        std::cout << "@ sizeof(int)  = " << sizeof(int)  << '\n';
        std::cout << "@ sizeof(int*) = " << sizeof(int*) << '\n';
        std::cout << "@ sizeof(long) = " << sizeof(long) << '\n';
        std::cout << "@ ---" << std::endl;
        }
      
      };
    
    static arma_first_extra_debug_message arma_first_extra_debug_message_run;
    }

#endif


//! @}
