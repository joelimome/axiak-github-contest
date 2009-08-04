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


//! \addtogroup op_stddev
//! @{


//! \brief
//! For each row or for each column, find the standard deviation.
//! The result is stored in a dense matrix that has either one column or one row.
//! The dimension for which the standard deviations are found is set via the stddev() function.
template<typename eT>
inline
void
op_stddev::apply(Mat<eT>& out, const Mat<eT>& X, const u32 norm_type, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (X.n_elem == 0), "op_stddev::apply(): given matrix has no elements" );
  
  arma_debug_check( (norm_type > 1), "op_stddev::apply(): incorrect usage. norm_type must be 0 or 1");
  arma_debug_check( (dim > 1),       "op_stddev::apply(): incorrect usage. dim must be 0 or 1"      );
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 0");
    
    out.set_size(1, X.n_cols);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      out[col] = std::sqrt( op_var::direct_var( X.colptr(col), X.n_rows, norm_type ) );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 1");
    
    out.set_size(X.n_rows, 1);
    
    const eT norm_val = (norm_type == 0) ? eT(X.n_cols-1) : eT(X.n_cols);
    
    for(u32 row=0; row<X.n_rows; ++row)
      {
      eT acc1 = eT(0);
      eT acc2 = eT(0);
  
      for(u32 col=0; col<X.n_cols; ++col)
        {
        const eT tmp_val = X.at(row,col);
        acc1 += tmp_val;
        acc2 += tmp_val*tmp_val;
        }
      
      const eT sd_val = std::sqrt( (acc2 - acc1*acc1/eT(X.n_cols)) / norm_val );
      
      out[row] = sd_val;
      }
    
    }
  
  }



//! implementation for complex numbers
template<typename T>
inline
void
op_stddev::apply(Mat<T>& out, const Mat< std::complex<T> >& X, const u32 norm_type, const u32 dim)
  {
  arma_extra_debug_sigprint();
  
  typedef typename std::complex<T> eT;
  
  arma_debug_check( (X.n_elem == 0), "op_stddev::apply(): given matrix has no elements" );
  
  arma_debug_check( (norm_type > 1), "op_stddev::apply(): incorrect usage. norm_type must be 0 or 1");
  arma_debug_check( (dim > 1),       "op_stddev::apply(): incorrect usage. dim must be 0 or 1"      );
  
  
  if(dim == 0)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 0");
    
    out.set_size(1, X.n_cols);
    
    for(u32 col=0; col<X.n_cols; ++col)
      {
      out[col] = std::sqrt( op_var::direct_var( X.colptr(col), X.n_rows, norm_type ) );
      }
    }
  else
  if(dim == 1)
    {
    arma_extra_debug_print("op_stddev::apply(), dim = 1");
    
    out.set_size(X.n_rows, 1);
    
    const T norm_val = (norm_type == 0) ? T(X.n_cols-1) : T(X.n_cols);
    
    for(u32 row=0; row<X.n_rows; ++row)
      {
      eT acc1 = eT(0);
      T  acc2 = T(0);
  
      for(u32 col=0; col<X.n_cols; ++col)
        {
        acc1 += X.at(row,col);
        acc2 += std::norm(X.at(row,col));
        }
      
      const T var_val = (acc2 - std::norm(acc1)/T(X.n_cols)) / norm_val;
      
      out[row] = std::sqrt(var_val);
      }
    
    }
  
  }



//! @}
