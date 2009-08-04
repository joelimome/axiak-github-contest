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

class arma_ostream_state
  {
  private:

  const ios::fmtflags   orig_flags;
  const std::streamsize orig_precision;
  const std::streamsize orig_width;
  const char            orig_fill;


  public:

  inline
  arma_ostream_state(const std::ostream& o)
   : orig_flags    (o.flags())
   , orig_precision(o.precision())
   , orig_width    (o.width())
   , orig_fill     (o.fill())
   {
   }


  inline
  void
  restore(std::ostream& o) const
    {
    o.flags    (orig_flags);
    o.precision(orig_precision);
    o.width    (orig_width);
    o.fill     (orig_fill);
    }

  };



class arma_ostream
  {
  public:
  
  template<typename eT>
  inline static u32 modify_stream(std::ostream& o, const Mat<eT>& m);
  
  template<typename T>
  inline static u32 modify_stream(std::ostream& o, const Mat< std::complex<T> >& m);
  
  template<typename eT>
  inline static void print(std::ostream& o, const Mat<eT>& m, const bool modify);

  template<typename T>
  inline static void print(std::ostream& o, const Mat< std::complex<T> >& m, const bool modify);
  };



template<typename eT>
inline
u32
arma_ostream::modify_stream(std::ostream& o, const Mat<eT>& m)
  {
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.fill(' ');

  u32 cell_width;
  
  bool use_layout_B = false;
  bool use_layout_C = false;
  
  for(u32 i=0; i<m.n_elem; ++i)
    {
    const eT val = m.mem[i];
    
    if(
      val >= eT(+100) ||
      ( (is_signed<eT>::value == true) && (val <= eT(-100)) ) ||
      ( (is_non_integral<eT>::value == true) && (val > eT(0)) && (val <= eT(+1e-4)) ) ||
      ( (is_non_integral<eT>::value == true) && (is_signed<eT>::value == true) && (val < eT(0)) && (val >= eT(-1e-4)) ) 
      )
      {
      use_layout_C = true;
      break;
      }
      
    if(
      (val >= eT(+10)) || ( (is_signed<eT>::value == true) && (val <= eT(-10)) )
      )
      {
      use_layout_B = true;
      }
    }
  
  if(use_layout_C == true)
    {
    o.setf(ios::scientific);
    o.unsetf(ios::fixed);
    o.precision(4);
    cell_width = 13;
    }
  else
  if(use_layout_B == true)
    {
    o.unsetf(ios::scientific);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 10;
    }
  else
    {
    o.unsetf(ios::scientific);
    o.setf(ios::fixed);
    o.precision(4);
    cell_width = 9;
    }
  
  return cell_width;
  }



//! "better than nothing" settings for complex numbers
template<typename T>
inline
u32
arma_ostream::modify_stream(std::ostream& o, const Mat< std::complex<T> >& m)
  {
  o.unsetf(ios::showbase);
  o.unsetf(ios::uppercase);
  o.fill(' ');
  
  o.setf(ios::scientific);
  o.setf(ios::showpos);
  o.unsetf(ios::fixed);
  
  u32 cell_width;
  
  o.precision(3);
  cell_width = 2 + 2*(1 + 3 + o.precision() + 5) + 1;

  return cell_width;
  }



//! Print a matrix to the specified stream
template<typename eT>
inline
void
arma_ostream::print(std::ostream& o, const Mat<eT>& m, const bool modify)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);

  u32 cell_width;
  
  if(modify == true)
    {
    cell_width = arma_ostream::modify_stream(o, m);
    }
  else
    {
    cell_width = o.width();
    }
  
  if(cell_width != 0)
    {
    for(u32 row=0; row < m.n_rows; ++row)
      {
      for(u32 col=0; col < m.n_cols; ++col)
        {
        o.width(cell_width);
        o << m.at(row,col);
        }
      
      o << '\n';
      }
    }
  else
    {
    for(u32 row=0; row < m.n_rows; ++row)
      {
      for(u32 col=0; col < m.n_cols-1; ++col)
        {
        o << m.at(row,col) << ' ';
        }
      
      o << m.at(row, m.n_cols-1) << '\n';
      }
    }
  
  o.flush();
  stream_state.restore(o);
  }



//! Print a complex matrix to the specified stream
//! EXPERIMENTAL !
template<typename T>
inline
void
arma_ostream::print(std::ostream& o, const Mat< std::complex<T> >& m, const bool modify)
  {
  arma_extra_debug_sigprint();
  
  const arma_ostream_state stream_state(o);

  u32 cell_width;
  
  if(modify == true)
    {
    cell_width = arma_ostream::modify_stream(o, m);
    }
  else
    {
    cell_width = o.width();
    }


  
  if(cell_width != 0)
    {
    for(u32 row=0; row < m.n_rows; ++row)
      {
      for(u32 col=0; col < m.n_cols; ++col)
        {
        std::ostringstream ss;
        ss.flags(o.flags());
        //ss.imbue(o.getloc());
        ss.precision(o.precision());

        ss << '(' << m.at(row,col).real() << ',' << m.at(row,col).imag() << ')';

        o.width(cell_width);
        o << ss.str();
        }

      o << '\n';
      }
    }
  else
    {
    for(u32 row=0; row < m.n_rows; ++row)
      {
      for(u32 col=0; col < m.n_cols-1; ++col)
        {
        o << '(' << m.at(row,col).real() << ',' << m.at(row,col).imag() << ") ";
        }
      o << '(' << m.at(row, m.n_cols-1).real() << ',' << m.at(row, m.n_cols-1).imag() << ")\n";
      }
    }

  
  o.flush();
  stream_state.restore(o);
  }



template<typename eT>
inline
std::ostream&
operator<< (std::ostream& o, const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  arma_ostream::print(o,m,true);
  
  return o;
  }



//! @}
