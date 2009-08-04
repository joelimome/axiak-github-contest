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


//! \addtogroup ostream
//! @{


//! Print a diagonal matrix to the specified stream.
template<typename T1>
inline
std::ostream&
operator<< (std::ostream& o, const Op<T1,op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1> tmp(X.m);
  const Mat<eT>& m = tmp.M;
    
  arma_debug_check( ((m.is_vec() == false) && (m.is_square() == false)), "operator<<(): incompatible dimensions for diagmat operation" );
  
  const arma_ostream_state stream_state(o);

  const u32 cell_width = arma_ostream::modify_stream(o, m);
  
  const u32 local_n_rows = (std::max)(m.n_rows, m.n_cols);
  
  for(u32 row=0; row < local_n_rows; ++row)
    {
    for(u32 col=0; col < local_n_rows; ++col)
      {
      if(row != col)
        {
        o.width(cell_width);
        if(is_complex<eT>::value == false)
          {
          o << "0.0";
          }
        else
          {
          o << "(0.0,0.0)";
          }
        }
      else
        {
        const eT val = m.is_vec() ? m.mem[row] : m.at(row,row);
        
        o.width(cell_width);
        o << val;
        }
      }
      o << '\n';
    }
  
  o.flush();

  stream_state.restore(o);  

  return o;
  }

//! @}
