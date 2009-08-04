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


//! \addtogroup static_assert
//! @{


//! Classes for primitive compile time assertions (until C++0x)
template<bool>
struct arma_static_assert;

template<>
struct arma_static_assert<true>
  {
  };


template<bool val>
struct arma_type_check
  {
  arma_inline static void apply()
    {
    arma_static_assert<!val> ERROR___INCORRECT_TYPE;
    ERROR___INCORRECT_TYPE = ERROR___INCORRECT_TYPE;
    }
  };


//
//
//

template<bool, class arma_class>
struct arma_apply_proxy;

template<class arma_class>
struct arma_apply_proxy<false, arma_class>
  {
  public:
  //template<typename T1> inline static void apply(const T1&);
  
  template<typename Tout, typename T1, typename T2>
  inline static void apply(Tout& out, const T1& A, const T2& B)
    {
    arma_static_assert<false> ERROR___INCORRECT_TYPE;
    ERROR___INCORRECT_TYPE = ERROR___INCORRECT_TYPE;
    }
  
  
  };

template<class arma_class>
struct arma_apply_proxy<true, arma_class>
  {
  template<typename Tout, typename T1, typename T2>
  inline static void apply(Tout& out, const T1& A, const T2& B)
    {
    arma_class::apply(out, A,B);
    }
  
  };


//! @}
