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


using std::cout;
using std::cerr;
using std::endl;
using std::ios;

template<typename eT> class Mat;
template<typename eT> class Col;
template<typename eT> class Row;

template<typename eT> class subview;
template<typename eT> class subview_col;
template<typename eT> class subview_row;
template<typename oT> class subview_field;

template<typename eT> class diagview;

class diskio;

class op_min;
class op_max;

class op_trans;
class op_htrans;
class op_conj;
class op_diagmat;
class op_inv;
class op_sum;
class op_neg;
class op_scalar_plus;
class op_scalar_minus_pre;
class op_scalar_minus_post;
class op_scalar_times;
class op_scalar_divide;

class glue_div;
class glue_minus;
class glue_plus;
class glue_times;
class glue_times_vec;
class glue_schur;

class glue_plus_diag;
class glue_minus_diag;
class glue_times_diag;
class glue_schur_diag;

template<const bool, const bool, const bool, const bool> class gemm;
template<const bool, const bool, const bool>       class gemv;

template<typename T1, typename op_type> class Op; 
template<typename T1, typename T2, typename glue_type> class Glue;



//! \addtogroup diskio
//! @{


//! file types supported by Armadillo
enum file_type
  {
  auto_detect,  //!< Automatically detect the file type (file must be one of the following types)
  raw_ascii,    //!< ASCII format (text), without any other information.
  arma_ascii,   //!< Armadillo ASCII format (text), with information about matrix type and size
  arma_binary,  //!< Armadillo binary format
  pgm_binary,   //!< Portable Grey Map (greyscale image)
  ppm_binary    //!< Portable Pixel Map (colour image), used by the field class only
  };


//! @}


