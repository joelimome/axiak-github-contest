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


//! \addtogroup Mat
//! @{


template<typename eT>
inline
Mat<eT>::~Mat()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(n_elem > sizeof(mem_local)/sizeof(eT) )
    {
    delete [] mem;
    }
    
  if(arma_config::debug == true)
    {
    // try to expose buggy user code that accesses deleted objects
    access::rw(n_rows) = 0;
    access::rw(n_cols) = 0;
    access::rw(n_elem) = 0;
    access::rw(mem)    = 0;
    }
  
  isnt_supported_elem_type<eT>::check();
  }



template<typename eT>
inline
Mat<eT>::Mat()
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  }


//! construct the matrix to have user specified dimensions
template<typename eT>
inline
Mat<eT>::Mat(const u32 in_n_rows, const u32 in_n_cols)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(in_n_rows, in_n_cols);
  }



//! change the matrix to have user specified dimensions (data is not preserved)
template<typename eT>
inline
void
Mat<eT>::set_size(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint(arma_boost::format("in_n_rows = %d, in_n_cols = %d") % in_n_rows % in_n_cols);
  
  init(in_n_rows,in_n_cols);
  }


//! internal matrix construction; if the requested size is small enough, memory from the stack is used. otherwise memory is allocated via 'new'
template<typename eT>
inline
void
Mat<eT>::init(const u32 in_n_rows, const u32 in_n_cols)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_n_rows = %d, in_n_cols = %d") % in_n_rows % in_n_cols );
  
  
  const u32 new_n_elem = in_n_rows * in_n_cols;

  if(n_elem == new_n_elem)
    {
    access::rw(n_rows) = in_n_rows;
    access::rw(n_cols) = in_n_cols;
    }
  else
    {
    
    if(n_elem > sizeof(mem_local)/sizeof(eT) )
      {
      delete [] mem;
      }
    
    if(new_n_elem <= sizeof(mem_local)/sizeof(eT) )
      {
      access::rw(mem) = mem_local;
      }
    else
      {
      access::rw(mem) = new(std::nothrow) eT[new_n_elem];
      arma_check( (mem == 0), "Mat::init(): out of memory" );
      }
    
    access::rw(n_elem) = new_n_elem;

    if(new_n_elem == 0)
      {
      access::rw(n_rows) = 0;
      access::rw(n_cols) = 0;
      }
    else
      {
      access::rw(n_rows) = in_n_rows;
      access::rw(n_cols) = in_n_cols;
      }
    
    }
  
  }


//! create the matrix from a textual description
template<typename eT>
inline
Mat<eT>::Mat(const char* text)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init( std::string(text) );
  }
  
  
  
//! create the matrix from a textual description
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const char* text)
  {
  arma_extra_debug_sigprint();
  
  init( std::string(text) );
  return *this;
  }
  
  

//! create the matrix from a textual description
template<typename eT>
inline
Mat<eT>::Mat(const std::string& text)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(text);
  }
  
  
  
//! create the matrix from a textual description
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  init(text);
  return *this;
  }



//! internal function to create the matrix from a textual description
template<typename eT>
inline 
void
Mat<eT>::init(const std::string& text)
  {
  arma_extra_debug_sigprint();
  
  //
  // work out the size
  
  u32 t_n_rows = 0;
  u32 t_n_cols = 0;
  
  bool t_n_cols_found = false;
  
  std::string token;
  
  std::string::size_type line_start = 0;
  std::string::size_type   line_end = 0;
  
  while( line_start < text.length() )
    {
    
    line_end = text.find(';', line_start);
    
    if(line_end == std::string::npos)
      line_end = text.length()-1;
    
    std::string::size_type line_len = line_end - line_start + 1;
    std::stringstream line_stream( text.substr(line_start,line_len) );
    
    
    u32 line_n_cols = 0;
    while(line_stream >> token)
      {
      ++line_n_cols;
      }
    
    
    if(line_n_cols > 0)
      {
      if(t_n_cols_found == false)
        {
        t_n_cols = line_n_cols;
        t_n_cols_found = true;
        }
      else
        arma_check( (line_n_cols != t_n_cols), "Mat::init(): inconsistent number of columns in given string");
      
      ++t_n_rows;
      }
    line_start = line_end+1;
    
    }
    
  Mat<eT> &x = *this;
  x.set_size(t_n_rows, t_n_cols);
  
  line_start = 0;
  line_end = 0;
  
  u32 row = 0;
  
  while( line_start < text.length() )
    {
    
    line_end = text.find(';', line_start);
    
    if(line_end == std::string::npos)
      line_end = text.length()-1;
    
    std::string::size_type line_len = line_end - line_start + 1;
    std::stringstream line_stream( text.substr(line_start,line_len) );
    
//     u32 col = 0;
//     while(line_stream >> token)
//       {
//       x.at(row,col) = strtod(token.c_str(), 0);
//       ++col;
//       }
    
    u32 col = 0;
    eT val;
    while(line_stream >> val)
      {
      x.at(row,col) = val;
      ++col;
      }
    
    ++row;
    line_start = line_end+1;
    }
  
  }



//! Set the matrix to be equal to the specified scalar.
//! NOTE: the size of the matrix will be 1x1
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  init(1,1);
  access::rw(mem[0]) = val;
  return *this;
  }



//! In-place addition of a scalar to all elements of the matrix
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) += val;
    }
  
  return *this;
  }



//! In-place subtraction of a scalar from all elements of the matrix
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) -= val;
    }
      
  return *this;
  }



//! In-place multiplication of all elements of the matrix with a scalar
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) *= val;
    }
  
  return *this;
  }



//! In-place division of all elements of the matrix with a scalar
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) /= val;
    }
  
  return *this;
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
Mat<eT>::Mat(const Mat<eT> &in_mat)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint(arma_boost::format("this = %x   in_mat = %x") % this % &in_mat);
  
  init(in_mat);
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const Mat<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  init(x);
  return *this;
  }



//! construct a matrix from a given matrix
template<typename eT>
inline
void
Mat<eT>::init(const Mat<eT> &x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    init(x.n_rows, x.n_cols);
    syslib::copy_elem( memptr(), x.mem, n_elem );
    }
  }



//! construct a matrix from a given auxillary array of eTs
template<typename eT>
inline
Mat<eT>::Mat(const eT* aux_mem, const u32 aux_n_rows, const u32 aux_n_cols)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(aux_n_rows, aux_n_cols);
  syslib::copy_elem( memptr(), aux_mem, n_elem );
  }



//! in-place matrix addition
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  glue_plus::apply_inplace(*this, m);
  return *this;
  }



//! in-place matrix subtraction
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  glue_minus::apply_inplace(*this, m);
  return *this;
  }



//! in-place matrix multiplication
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator*=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  glue_times::apply_inplace(*this, m);
  return *this;
  }



//! in-place element-wise matrix multiplication
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator%=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  glue_schur::apply_inplace(*this, m);
  return *this;
  }



//! in-place element-wise matrix division
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator/=(const Mat<eT>& m)
  {
  arma_extra_debug_sigprint();
  
  glue_div::apply_inplace(*this, m);
  return *this;
  }



//! for constructing a complex matrix out of two non-complex matrices
template<typename eT>
template<typename T1, typename T2>
inline
Mat<eT>::Mat
  (
  const Base<typename Mat<eT>::pod_type,T1>& A,
  const Base<typename Mat<eT>::pod_type,T2>& B
  )
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  arma_type_check< is_complex<eT>::value == false >::apply();   //!< compile-time abort if eT isn't std::complex
  
  typedef typename T1::elem_type T;
  arma_type_check< is_complex<T>::value == true >::apply();   //!< compile-time abort if T is std::complex
  
  isnt_same_type<std::complex<T>, eT>::check();   //!< compile-time abort if types are not compatible
  
  const unwrap<T1> tmp_A(A.get_ref());
  const unwrap<T2> tmp_B(B.get_ref());
  
  const Mat<T>& X = tmp_A.M;
  const Mat<T>& Y = tmp_B.M;
  
  arma_assert_same_size(X, Y, "Mat()");
  
  init(X.n_rows, Y.n_cols);
  
  const T* X_mem = X.mem;
  const T* Y_mem = Y.mem;
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) = std::complex<T>(X_mem[i], Y_mem[i]);
    }
  }



//! construct a matrix from subview (e.g. construct a matrix from a delayed submatrix operation)
template<typename eT>
inline
Mat<eT>::Mat(const subview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a matrix from subview (e.g. construct a matrix from a delayed submatrix operation)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::extract(*this, X);
  return *this;
  }


//! in-place matrix addition (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator+=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::plus_inplace(*this, X);
  return *this;
  }


//! in-place matrix subtraction (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator-=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::minus_inplace(*this, X);
  return *this;
  }



//! in-place matrix mutiplication (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator*=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::times_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix mutiplication (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator%=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::schur_inplace(*this, X);
  return *this;
  }



//! in-place element-wise matrix division (using a submatrix on the right-hand-side)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator/=(const subview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  subview<eT>::div_inplace(*this, X);
  return *this;
  }



//! construct a matrix from diagview (e.g. construct a matrix from a delayed diag operation)
template<typename eT>
inline
Mat<eT>::Mat(const diagview<eT>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  
  this->operator=(X);
  }



//! construct a matrix from diagview (e.g. construct a matrix from a delayed diag operation)
template<typename eT>
inline
const Mat<eT>&
Mat<eT>::operator=(const diagview<eT>& X)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>::extract(*this, X);
  return *this;
  }



//! creation of subview (row vector)
template<typename eT>
arma_inline
subview_row<eT>
Mat<eT>::row(const u32 row_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "Mat::row(): row out of bounds" );
  
  return subview_row<eT>(*this, row_num);
  }



//! creation of subview (row vector)
template<typename eT>
arma_inline
const subview_row<eT>
Mat<eT>::row(const u32 row_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( row_num >= n_rows, "Mat::row(): row out of bounds" );
  
  return subview_row<eT>(*this, row_num);
  }



//! creation of subview (column vector)
template<typename eT>
arma_inline
subview_col<eT>
Mat<eT>::col(const u32 col_num)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "Mat::col(): out of bounds");
  
  return subview_col<eT>(*this, col_num);
  }



//! creation of subview (column vector)
template<typename eT>
arma_inline
const subview_col<eT>
Mat<eT>::col(const u32 col_num) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( col_num >= n_cols, "Mat::col(): out of bounds");
  
  return subview_col<eT>(*this, col_num);
  }



//! creation of subview (submatrix comprised of specified row vectors)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "Mat::rows(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, in_row1, 0, in_row2, n_cols-1);
  }



//! creation of subview (submatrix comprised of specified row vectors)
template<typename eT>
arma_inline
const subview<eT>
Mat<eT>::rows(const u32 in_row1, const u32 in_row2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_row2 >= n_rows),
    "Mat::rows(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, in_row1, 0, in_row2, n_cols-1);
  }



//! creation of subview (submatrix comprised of specified column vectors)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "Mat::cols(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, 0, in_col1, n_rows-1, in_col2);
  }



//! creation of subview (submatrix comprised of specified column vectors)
template<typename eT>
arma_inline
const subview<eT>
Mat<eT>::cols(const u32 in_col1, const u32 in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 > in_col2) || (in_col2 >= n_cols),
    "Mat::cols(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, 0, in_col1, n_rows-1, in_col2);
  }



//! creation of subview (submatrix)
template<typename eT>
arma_inline
subview<eT>
Mat<eT>::submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "Mat::submat(): indices out of bounds or incorrectly used"
    );
  
  return subview<eT>(*this, in_row1, in_col1, in_row2, in_col2);
  }



//! creation of subview (generic submatrix)
template<typename eT>
arma_inline
const subview<eT>
Mat<eT>::submat(const u32 in_row1, const u32 in_col1, const u32 in_row2, const u32 in_col2) const
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 > in_row2) || (in_col1 >  in_col2) || (in_row2 >= n_rows) || (in_col2 >= n_cols),
    "Mat::submat(): indices out of bounds or incorrectly used"
    );
    
  return subview<eT>(*this, in_row1, in_col1, in_row2, in_col2);
  }



//! creation of diagview (diagonal)
template<typename eT>
arma_inline
diagview<eT>
Mat<eT>::diag(const s32 in_id)
  {
  arma_extra_debug_sigprint();
  
  const u32 row_offset = (in_id < 0) ? -in_id : 0;
  const u32 col_offset = (in_id > 0) ?  in_id : 0;
  
  arma_debug_check
    (
    (row_offset >= n_rows) || (col_offset >= n_cols),
    "Mat::diag(): out of bounds"
    );
  
  const u32 len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  return diagview<eT>(*this, row_offset, col_offset, len);
  }



//! creation of diagview (diagonal)
template<typename eT>
arma_inline
const diagview<eT>
Mat<eT>::diag(const s32 in_id) const
  {
  arma_extra_debug_sigprint();
  
  const u32 row_offset = (in_id < 0) ? -in_id : 0;
  const u32 col_offset = (in_id > 0) ?  in_id : 0;
  
  arma_debug_check
    (
    (row_offset >= n_rows) || (col_offset >= n_cols),
    "Mat::diag(): out of bounds"
    );
  
  
  const u32 len = (std::min)(n_rows - row_offset, n_cols - col_offset);
  
  return diagview<eT>(*this, row_offset, col_offset, len);
  }



template<typename eT>
inline
void
Mat<eT>::swap_rows(const u32 in_row1, const u32 in_row2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_row1 >= n_rows) || (in_row2 >= n_rows),
    "Mat::swap_rows(): out of bounds"
    );
  
  for(u32 col=0; col<n_cols; ++col)
    {
    const u32 offset = col*n_rows;
    const u32 pos1   = in_row1 + offset;
    const u32 pos2   = in_row2 + offset;
    
    const eT tmp          = mem[pos1];
    access::rw(mem[pos1]) = mem[pos2];
    access::rw(mem[pos2]) = tmp;
    }
  
  }



template<typename eT>
inline
void
Mat<eT>::swap_cols(const u32 in_col1, const u32 in_col2)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check
    (
    (in_col1 >= n_cols) || (in_col2 >= n_cols),
    "Mat::swap_cols(): out of bounds"
    );
  
  eT* ptr1 = colptr(in_col1);
  eT* ptr2 = colptr(in_col2);
  
  for(u32 row=0; row<n_rows; ++row)
    {
    const eT tmp = ptr1[row];
    ptr1[row]    = ptr2[row];
    ptr2[row]    = tmp;
    }
  
  }



//! create a matrix from Op, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename op_type>
inline
Mat<eT>::Mat(const Op<T1, op_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  op_type::apply(*this, X);
  }



//! create a matrix from Op, i.e. run the previously delayed unary operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, typename T1::elem_type>::check();
  
  op_type::apply(*this, X);
  
  return *this;
  }



//! in-place matrix addition, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  glue_plus::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix subtraction, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  glue_minus::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix multiplication, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  
  glue_times::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix element-wise multiplication, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  glue_schur::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix element-wise division, with the right-hand-side operand having delayed operations
template<typename eT>
template<typename T1, typename op_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const Op<T1, op_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  glue_div::apply_inplace(*this, X);
  
  return *this;
  }



//! create a matrix from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
Mat<eT>::Mat(const Glue<T1, T2, glue_type>& X)
  : n_rows(0)
  , n_cols(0)
  , n_elem(0)
  //, mem(0)
  , mem(mem)
  {
  arma_extra_debug_sigprint_this(this);
  this->operator=(X);
  }



//! create a matrix from Glue, i.e. run the previously delayed binary operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  // TODO:
  // it may be simpler to pass the two objects (currently wrapped in X)
  // directly to the apply function.
  // (many adjustments throughout the source code will be required)
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_type::apply(*this, X);
  
  return *this;
  }


//! in-place matrix addition, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator+=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_plus::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix subtraction, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator-=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_minus::apply_inplace(*this, X);
  
  return *this;
  }



//! in-place matrix multiplications, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator*=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_times::apply_inplace(*this, X);
  return *this;
  }



//! in-place matrix element-wise multiplication, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator%=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_schur::apply_inplace(*this, X);
  return *this;
  }



//! in-place matrix element-wise division, with the right-hand-side operands having delayed operations
template<typename eT>
template<typename T1, typename T2, typename glue_type>
inline
const Mat<eT>&
Mat<eT>::operator/=(const Glue<T1, T2, glue_type>& X)
  {
  arma_extra_debug_sigprint();
  
  isnt_same_type<eT, typename T1::elem_type>::check();
  isnt_same_type<eT, typename T2::elem_type>::check();
  
  glue_div::apply_inplace(*this, X);
  return *this;
  }



//! linear element accessor (treats the matrix as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT&
Mat<eT>::operator() (const u32 i)
  {
  arma_debug_check( (i >= n_elem), "Mat::operator(): index out of bounds");
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the matrix as a vector); bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT
Mat<eT>::operator() (const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "Mat::operator(): index out of bounds");
  return mem[i];
  }


//! linear element accessor (treats the matrix as a vector); no bounds check.  
template<typename eT>
arma_inline
eT&
Mat<eT>::operator[] (const u32 i)
  {
  return access::rw(mem[i]);
  }



//! linear element accessor (treats the matrix as a vector); no bounds check
template<typename eT>
arma_inline
eT
Mat<eT>::operator[] (const u32 i) const
  {
  return mem[i];
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT&
Mat<eT>::operator() (const u32 in_row, const u32 in_col)
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "Mat::operator(): index out of bounds");
  return access::rw(mem[in_row + in_col*n_rows]);
  }



//! element accessor; bounds checking not done when ARMA_NO_DEBUG is defined
template<typename eT>
arma_inline
eT
Mat<eT>::operator() (const u32 in_row, const u32 in_col) const
  {
  arma_debug_check( ((in_row >= n_rows) || (in_col >= n_cols)), "Mat::operator(): index out of bounds");
  return mem[in_row + in_col*n_rows];
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
eT&
Mat<eT>::at(const u32 in_row, const u32 in_col)
  {
  return access::rw( mem[in_row + in_col*n_rows] );
  }



//! element accessor; no bounds check
template<typename eT>
arma_inline
eT
Mat<eT>::at(const u32 in_row, const u32 in_col) const
  {
  return mem[in_row + in_col*n_rows];
  }



//! prefix ++
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator++()
  {
  Mat_aux::prefix_pp(*this);
  return *this;
  }



//! postfix ++  (must not return the object by reference)
template<typename eT>
arma_inline
void
Mat<eT>::operator++(int)
  {
  Mat_aux::postfix_pp(*this);
  }



//! prefix --
template<typename eT>
arma_inline
const Mat<eT>&
Mat<eT>::operator--()
  {
  Mat_aux::prefix_mm(*this);
  return *this;
  }



//! postfix --  (must not return the object by reference)
template<typename eT>
arma_inline
void
Mat<eT>::operator--(int)
  {
  Mat_aux::postfix_mm(*this);
  }



//! returns true if the object can be interpreted as a column or row vector
template<typename eT>
arma_inline
bool
Mat<eT>::is_vec() const
  {
  return ( (n_rows == 1) || (n_cols == 1) );
  }



//! returns true if the object has the same number of non-zero rows and columnns
template<typename eT>
arma_inline
bool
Mat<eT>::is_square() const
  {
  return ( (n_rows == n_cols) && (n_elem > 0) );
  }



//! returns true if all of the elements are finite
template<typename eT>
arma_inline
bool
Mat<eT>::is_finite() const
  {
  for(u32 i=0; i<n_elem; ++i)
    {
    if(arma_isfinite(mem[i]) == false)
      {
      return false;
      }
    }

  return true;
  }



//! returns a pointer to array of eTs for a specified column; no bounds check
template<typename eT>
arma_inline
eT*
Mat<eT>::colptr(const u32 in_col)
  {
  return & access::rw(mem[in_col*n_rows]);
  }



//! returns a pointer to array of eTs for a specified column; no bounds check
template<typename eT>
arma_inline
const eT*
Mat<eT>::colptr(const u32 in_col) const
  {
  return & mem[in_col*n_rows];
  }



//! returns a pointer to array of eTs used by the matrix
template<typename eT>
arma_inline
eT*
Mat<eT>::memptr()
  {
  return const_cast<eT*>(mem);
  }



//! returns a pointer to array of eTs used by the matrix
template<typename eT>
arma_inline
const eT*
Mat<eT>::memptr() const
  {
  return mem;
  }



//! print contents of the matrix (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Mat<eT>::print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    cout << extra_text << '\n';
    }
  
  arma_ostream::print(cout, *this, true);
  }


//! print contents of the matrix to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the precision and cell width are modified.
//! on return, the stream's flags are restored to their original values.
template<typename eT>
inline
void
Mat<eT>::print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    user_stream << extra_text << '\n';
    }
  
  arma_ostream::print(user_stream, *this, true);
  }



//! print contents of the matrix (to the cout stream),
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Mat<eT>::raw_print(const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    cout << extra_text << '\n';
    }
  
  arma_ostream::print(cout, *this, false);
  }



//! print contents of the matrix to a user specified stream,
//! optionally preceding with a user specified line of text.
//! the stream's flags are used as is and are not modified.
//! (i.e. the precision and cell width are not modified).
template<typename eT>
inline
void
Mat<eT>::raw_print(std::ostream& user_stream, const std::string extra_text) const
  {
  arma_extra_debug_sigprint();
  
  if(extra_text.length() != 0)
    {
    user_stream << extra_text << '\n';
    }
  
  arma_ostream::print(user_stream, *this, false);
  }



//! fill the matrix with the specified value
template<typename eT>
inline
void
Mat<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) = val;
    }
  }



template<typename eT>
inline
void
Mat<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  fill(eT(0));
  }



template<typename eT>
inline
void
Mat<eT>::zeros(const u32 in_rows, const u32 in_cols)
  {
  arma_extra_debug_sigprint( arma_boost::format("in_rows = %d, in_cols = %d") % in_rows % in_cols );

  set_size(in_rows, in_cols);
  fill(eT(0));
  }



template<typename eT>
inline
void
Mat<eT>::reset()
  {
  arma_extra_debug_sigprint();
  
  init(0,0);
  }



//! save the matrix to a file
template<typename eT>
inline
void
Mat<eT>::save(const std::string name, const file_type type) const
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case raw_ascii:
      diskio::save_raw_ascii(*this, name);
      break;
    
    case arma_ascii:
      diskio::save_arma_ascii(*this, name);
      break;
    
    case arma_binary:
      diskio::save_arma_binary(*this, name);
      break;
      
    case pgm_binary:
      diskio::save_pgm_binary(*this, name);
      break;
    
    default:
      arma_stop("Mat::save(): unsupported type");
    }
  
  }



//! load a matrix from a file
template<typename eT>
inline
void
Mat<eT>::load(const std::string name, const file_type type)
  {
  arma_extra_debug_sigprint();
  
  switch(type)
    {
    case auto_detect:
      diskio::load_auto_detect(*this, name);
      break;
    
    case raw_ascii:
      diskio::load_raw_ascii(*this, name);
      break;
    
    case arma_ascii:
      diskio::load_arma_ascii(*this, name);
      break;
    
    case arma_binary:
      diskio::load_arma_binary(*this, name);
      break;
      
    case pgm_binary:
      diskio::load_pgm_binary(*this, name);
      break;
    
    default:
      arma_stop("Mat::load(): unsupported type");
    }
  
  }



//! prefix ++
template<typename eT>
arma_inline
void
Mat_aux::prefix_pp(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;

  for(u32 i=0; i<n_elem; ++i)
    {
    ++(memptr[i]);
    }
  }



//! prefix ++ for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::prefix_pp(Mat< std::complex<T> >& x)
  {
  x += T(1);
  }



//! postfix ++
template<typename eT>
arma_inline
void
Mat_aux::postfix_pp(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;

  for(u32 i=0; i<n_elem; ++i)
    {
    (memptr[i])++;
    }
  }



//! postfix ++ for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::postfix_pp(Mat< std::complex<T> >& x)
  {
  x += T(1);
  }



//! prefix --
template<typename eT>
arma_inline
void
Mat_aux::prefix_mm(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;

  for(u32 i=0; i<n_elem; ++i)
    {
    --(memptr[i]);
    }
  }



//! prefix -- for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::prefix_mm(Mat< std::complex<T> >& x)
  {
  x -= T(1);
  }



//! postfix --
template<typename eT>
arma_inline
void
Mat_aux::postfix_mm(Mat<eT>& x)
  {
        eT* memptr = x.memptr();
  const u32 n_elem = x.n_elem;

  for(u32 i=0; i<n_elem; ++i)
    {
    (memptr[i])--;
    }
  }



//! postfix ++ for complex numbers (work around for limitations of the std::complex class)
template<typename T>
arma_inline
void
Mat_aux::postfix_mm(Mat< std::complex<T> >& x)
  {
  x -= T(1);
  }



//! @}
