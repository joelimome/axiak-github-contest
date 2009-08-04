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


//! \addtogroup op_randn
//! @{


// TODO:
// implement a more efficient method
//
// possible option:
// Marsaglia and Tsang Ziggurat technique to transform from a uniform to a normal distribution.
// G. Marsaglia, W.W. Tsang.
// "Ziggurat method for generating random variables",
// J. Statistical Software, vol 5, 2000.
// http://www.jstatsoft.org/v05/i08/


template<typename eT>
inline
eT
op_randn::randn()
  {
  const u32 N = 12;
  eT acc = eT(0);
  
  for(u32 i=0; i<N; ++i)
    {
    acc += eT(std::rand()) / eT(RAND_MAX);
    }
  
  return acc/eT(N) - eT(0.5);
  }



template<typename eT>
inline
void
op_randn::direct_randn(eT* x, const u32 n_elem)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    x[i] = eT(op_randn::randn<eT>());
    }
  
  }



template<typename eT>
inline
void
op_randn::apply(Mat<eT>& out, const Op<Mat<eT>,op_randn>& in)
  {
  arma_extra_debug_sigprint();
  
  const u32 n_rows = in.aux_u32_a;
  const u32 n_cols = (in.aux_u32_b > 0) ? in.aux_u32_b : 1;
  
  out.set_size(n_rows, n_cols);
  
  op_randn::direct_randn(out.memptr(), out.n_elem);
  }



template<typename eT>
inline
void
op_randn::apply(Mat<eT>& out, const Op<Col<eT>,op_randn>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a, 1);
  
  op_randn::direct_randn(out.memptr(), out.n_elem);
  }



template<typename eT>
inline
void
op_randn::apply(Mat<eT>& out, const Op<Row<eT>,op_randn>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(1, in.aux_u32_a);
  
  op_randn::direct_randn(out.memptr(), out.n_elem);
  }



template<typename eT>
inline
void
op_randn::apply(Col<eT>& out, const Op<Col<eT>,op_randn>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a);
  
  op_randn::direct_randn(out.memptr(), out.n_elem);
  }



template<typename eT>
inline
void
op_randn::apply(Row<eT>& out, const Op<Row<eT>,op_randn>& in)
  {
  arma_extra_debug_sigprint();
  
  out.set_size(in.aux_u32_a);
  
  op_randn::direct_randn(out.memptr(), out.n_elem);
  }



//! @}
