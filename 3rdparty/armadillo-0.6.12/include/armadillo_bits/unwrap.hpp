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


//! \addtogroup unwrap
//! @{


template<typename T1>
class unwrap
  {
  public:
  inline unwrap(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



//template<>
template<typename eT>
class unwrap< Mat<eT> >
  {
  public:

  inline unwrap(const Mat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Mat<eT>& M;
  
  };



//template<>
template<typename eT>
class unwrap< Row<eT> >
  {
  public:

  inline unwrap(const Row<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Row<eT>& M;
  
  };



//template<>
template<typename eT>
class unwrap< Col<eT> >
  {
  public:

  inline unwrap(const Col<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Col<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap< Op<T1, op_type> >
  {
  public:

  inline unwrap(const Op<T1, op_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  typedef typename T1::elem_type elem_type;
  const Mat<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_type>
class unwrap< Glue<T1, T2, glue_type> >
  {
  public:

  inline unwrap(const Glue<T1, T2, glue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  typedef typename T1::elem_type elem_type;
  const Mat<elem_type> M;
  
  };


//template<>
template<typename eT>
class unwrap< subview<eT> >
  {
  public:

  inline unwrap(const subview<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Mat<eT> M;
  
  };


//template<>
template<typename eT>
class unwrap< diagview<eT> >
  {
  public:

  inline unwrap(const diagview<eT> & A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Mat<eT> M;
  
  };



//
//
//


template<typename T1>
class unwrap_to_elem_access
  {
  public:
  inline unwrap_to_elem_access(const T1& A)
    {
    arma_type_check< is_arma_type<T1>::value == false >::apply();
    }
  };



//template<>
template<typename eT>
class unwrap_to_elem_access< Mat<eT> >
  {
  public:

  inline unwrap_to_elem_access(const Mat<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  
  const Mat<eT>& M;
  
  };



//template<>
template<typename eT>
class unwrap_to_elem_access< Row<eT> >
  {
  public:

  inline unwrap_to_elem_access(const Row<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const Row<eT>& M;
  
  };



//template<>
template<typename eT>
class unwrap_to_elem_access< Op< Row<eT>, op_trans> >
  {
  public:

  // NOTE:
  // currently accessing M.n_rows and M.n_cols will give wrong results
  
  inline unwrap_to_elem_access(const Op<Row<eT>, op_trans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const Row<eT>& M;
  
  };



//template<>
template<typename eT>
class unwrap_to_elem_access< Col<eT> >
  {
  public:

  inline unwrap_to_elem_access(const Col<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const Col<eT>& M;
  
  };



//template<>
template<typename eT>
class unwrap_to_elem_access< Op<Col<eT>, op_trans> >
  {
  public:

  inline unwrap_to_elem_access(const Op<Col<eT>, op_trans>& A)
    : M(A.m)
    {
    arma_extra_debug_sigprint();
    }

  inline eT  operator[](const u32 i) const { return M[i]; }
  inline eT  operator()(const u32 i) const { return M(i); }
  
  // NOTE: use of in_row and in_col is swapped (due to transpose operation)
  inline eT          at(const u32 in_row, const u32 in_col) const { return M.at(in_col,in_row); }
  inline eT  operator()(const u32 in_row, const u32 in_col) const { return M(in_col,in_row);    }
  
  const Col<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap_to_elem_access< Op<T1, op_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline unwrap_to_elem_access(const Op<T1, op_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline elem_type  operator[](const u32 i) const { return M[i]; }
  inline elem_type  operator()(const u32 i) const { return M(i); }
  
  inline elem_type  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline elem_type          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const Mat<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_type>
class unwrap_to_elem_access< Glue<T1, T2, glue_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline unwrap_to_elem_access(const Glue<T1, T2, glue_type>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline elem_type  operator[](const u32 i) const { return M[i]; }
  inline elem_type  operator()(const u32 i) const { return M(i); }
  
  inline elem_type  operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline elem_type          at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
    
  const Mat<elem_type> M;
  
  };


//template<>
template<typename eT>
class unwrap_to_elem_access< subview<eT> >
  {
  public:

  inline unwrap_to_elem_access(const subview<eT>& A)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  inline eT operator[](const u32 i) const { return M[i]; }
  inline eT operator()(const u32 i) const { return M(i); }
  
  inline eT operator()(const u32 in_row, const u32 in_col) const { return M(in_row,in_col);    }
  inline eT         at(const u32 in_row, const u32 in_col) const { return M.at(in_row,in_col); }
  
  const subview<eT>& M;
  
  };


//
//
//

template<typename T1>
class unwrap_check
  {
  private:
  template<typename eT> inline unwrap_check(const T1& A, const Mat<eT>& B);
  template<typename eT> inline unwrap_check(const T1& A, const Row<eT>& B);
  template<typename eT> inline unwrap_check(const T1& A, const Col<eT>& B);
  template<typename eT> inline unwrap_check(const T1& A, const subview<eT>&  B);
  template<typename eT> inline unwrap_check(const T1& A, const diagview<eT>& B);
  };


//template <>
template<typename eT>
class unwrap_check< Mat<eT> >
  {
  public:

  inline
  unwrap_check(const Mat<eT>& A, const Mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const Mat<eT>*>(&B)) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Mat<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Mat<eT>& A, const Row<eT>& B)
    : M_local( (&A == reinterpret_cast<const Mat<eT>*>(&B)) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Mat<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Mat<eT>& A, const Col<eT>& B)
    : M_local( (&A == reinterpret_cast<const Mat<eT>*>(&B)) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Mat<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const Mat<eT>& A, const subview<eT>& B)
    : M_local( (&A == reinterpret_cast<const Mat<eT>*>(&B.m)) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Mat<eT>*>(&B.m)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const Mat<eT>& A, const diagview<eT>& B)
    : M_local( (&A == reinterpret_cast<const Mat<eT>*>(&B.m)) ? new Mat<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Mat<eT>*>(&B.m)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      {
      delete M_local;
      }
    }
  
  
  // the order below is important
  const Mat<eT>* M_local;
  const Mat<eT>& M;
  
  };



//template <>
template<typename eT>
class unwrap_check< Row<eT> >
  {
  public:

  inline
  unwrap_check(const Row<eT>& A, const Mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? new Row<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Row<eT>& A, const Row<eT>& B)
    : M_local( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? new Row<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  inline
  unwrap_check(const Row<eT>& A, const Col<eT>& B)
    : M_local( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? new Row<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Row<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const Row<eT>& A, const subview<eT>& B)
    : M_local( (&A == reinterpret_cast<const Row<eT>*>(&B.m)) ? new Row<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Row<eT>*>(&B.m)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const Row<eT>& A, const diagview<eT>& B)
    : M_local( (&A == reinterpret_cast<const Row<eT>*>(&B.m)) ? new Row<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Row<eT>*>(&B.m)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }

  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      delete M_local;
    }
  
  
  // the order below is important
  const Row<eT>* M_local;
  const Row<eT>& M;
  
  };




//template <>
template<typename eT>
class unwrap_check< Col<eT> >
  {
  public:

  inline
  unwrap_check(const Col<eT>& A, const Mat<eT>& B)
    : M_local( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? new Col<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Col<eT>& A, const Row<eT>& B)
    : M_local( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? new Col<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Col<eT>& A, const Col<eT>& B)
    : M_local( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? new Col<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Col<eT>*>(&B)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const Col<eT>& A, const subview<eT>& B)
    : M_local( (&A == reinterpret_cast<const Col<eT>*>(&B.m)) ? new Col<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Col<eT>*>(&B.m)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }


  inline
  unwrap_check(const Col<eT>& A, const diagview<eT>& B)
    : M_local( (&A == reinterpret_cast<const Col<eT>*>(&B.m)) ? new Col<eT>(A) : 0 )
    , M      ( (&A == reinterpret_cast<const Col<eT>*>(&B.m)) ? (*M_local)     : A )
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    if(M_local)
      delete M_local;
    }
  
  
  // the order below is important
  const Col<eT>* M_local;
  const Col<eT>& M;
  
  };



template<typename T1, typename op_type>
class unwrap_check< Op<T1, op_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  //template<typename eT>
  inline
  unwrap_check(const Op<T1,op_type>& A, const Mat<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  //template<typename eT>
  inline
  unwrap_check(const Op<T1,op_type>& A, const Row<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  //template<typename eT>
  inline
  unwrap_check(const Op<T1,op_type>& A, const Col<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    }
  
  const Mat<elem_type> M;
  
  };



template<typename T1, typename T2, typename glue_type>
class unwrap_check< Glue<T1, T2, glue_type> >
  {
  public:
  typedef typename T1::elem_type elem_type;

  inline
  unwrap_check(const Glue<T1, T2, glue_type>& A, const Mat<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Glue<T1, T2, glue_type>& A, const Row<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  unwrap_check(const Glue<T1, T2, glue_type>& A, const Col<elem_type>& B)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }
  
  
  inline
  ~unwrap_check()
    {
    arma_extra_debug_sigprint();
    
    }
  
  
  const Mat<elem_type> M;
  
  };




//template<>
template<typename eT>
class unwrap_check< subview<eT> >
  {
  public:

  template<typename T2>
  inline unwrap_check(const subview<eT>& A, const T2& junk)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Mat<eT> M;
  
  };


//template<>
template<typename eT>
class unwrap_check< diagview<eT> >
  {
  public:

  template<typename T2>
  inline unwrap_check(const diagview<eT> & A, const T2& junk)
    : M(A)
    {
    arma_extra_debug_sigprint();
    }

  const Mat<eT> M;
  
  };


//! @}
