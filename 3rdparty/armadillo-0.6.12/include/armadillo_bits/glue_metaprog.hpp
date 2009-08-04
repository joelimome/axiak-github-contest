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


//! \addtogroup glue_metaprog
//! @{


//! \brief
//! Template metaprogram depth_lhs
//! calculates the number of Glue<Tx,Ty, glue_type> instances on the left hand side argument of Glue<Tx,Ty, glue_type>
//! i.e. it recursively expands each Tx, until the type of Tx is not "Glue<..,.., glue_type>"  (i.e the "glue_type" changes)

template<typename glue_type, typename T1>
struct depth_lhs
  {
  static const u32 num = 0;
  };

template<typename glue_type, typename T1, typename T2>
struct depth_lhs< glue_type, Glue<T1,T2,glue_type> >
  {
  static const u32 num = 1 + depth_lhs<glue_type, T1>::num;
  };



//! \brief
//! Template metaprogram mat_ptrs
//! fills a given array with addresses of matrices from a recursive instance of Glue<Tx,Ty, glue_type>.
//! While parsing the recursive instance, if encountered objects are of type Op<..>,
//! they are converted to type 'Mat' first

template<typename glue_type, typename T1>
struct mat_ptrs
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 0;

  inline
  static
  void
  get_ptrs
    (
    const Mat<elem_type>** ptrs,
    bool* del,
    const T1& X
    )
    {

    ptrs[0] = 
      (
      is_Mat<T1>::value ?
        reinterpret_cast<const Mat<elem_type>*>(&X)
      :
        new Mat<elem_type>(X)
      );

    
    del[0] = 
      (
      is_Mat<T1>::value ?
        false
      :
        true
      );

    
    }
  
  };



template<typename glue_type, typename T1, typename T2>
struct mat_ptrs<glue_type, Glue<T1,T2,glue_type> >
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 1 + mat_ptrs<glue_type, T1>::num;
  
  inline
  static
  void
  get_ptrs
    (
    const Mat<elem_type>** in_ptrs,
    bool* del,
    const Glue<T1,T2,glue_type>& X
    )
    {
    isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
    
    mat_ptrs<glue_type, T1>::get_ptrs(in_ptrs, del, X.A);
    
    in_ptrs[num]  = 
      (
      is_Mat<T2>::value ?
        reinterpret_cast<const Mat<elem_type>*>(&X.B)
      :
        new Mat<elem_type>(X.B)
      );
    
    del[num] = 
      (
      is_Mat<T2>::value ?
        false
      :
        true
      );
    }
  
  };



//! template metaprogram mat_ptrs_outcheck
//! builds on 'mat_ptrs' by also checking whether any of the input matrices are aliases of the output matrix

template<typename glue_type, typename T1>
struct mat_ptrs_outcheck
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 0;

  inline
  static
  void
  get_ptrs
    (
    const Mat<elem_type>** ptrs,
    bool* del,
    const T1& X,
    const Mat<elem_type>* out_ptr
    )
    {

    const bool same_ptr = 
      (
      is_Mat<T1>::value ?
        (
        (out_ptr == reinterpret_cast<const Mat<elem_type>*>(&X)) ?
          true
        :
          false
        )
      :
        false
      );

    
    ptrs[0] = 
      (
      same_ptr ?
        new Mat<elem_type>(X)
      :
        (
        is_Mat<T1>::value ?
          reinterpret_cast<const Mat<elem_type>*>(&X)
        :
          new Mat<elem_type>(X)
        )
      );

    
    del[0] = 
      (
      same_ptr ?
        true
      :
        (
        is_Mat<T1>::value ?
          false
        :
          true
        )
      );

    
    }
  
  };



template<typename glue_type, typename T1, typename T2>
struct mat_ptrs_outcheck<glue_type, Glue<T1,T2,glue_type> >
  {
  typedef typename T1::elem_type elem_type;
  
  static const u32 num = 1 + mat_ptrs_outcheck<glue_type, T1>::num;
  
  inline
  static
  void
  get_ptrs
    (
    const Mat<elem_type>** ptrs,
    bool* del,
    const Glue<T1,T2,glue_type>& X,
    const Mat<elem_type>* out_ptr
    )
    {
    isnt_same_type<typename T1::elem_type, typename T2::elem_type>::check();
    
    mat_ptrs_outcheck<glue_type, T1>::get_ptrs(ptrs, del, X.A, out_ptr);
    
    const bool same_ptr =
      (
      is_Mat<T2>::value ?
        (
        (out_ptr == reinterpret_cast<const Mat<elem_type>*>(&X.B)) ?
          true
        :
          false
        )
      :
        false
      );
    
    
    ptrs[num]  = 
      (
      same_ptr ?
        new Mat<elem_type>(X.B)
      :
        (
        is_Mat<T2>::value ?
          reinterpret_cast<const Mat<elem_type>*>(&X.B)
        :
          new Mat<elem_type>(X.B)
        )
      );
    
    
    del[num] = 
      (
      same_ptr ?
        true
      :
        (
        is_Mat<T2>::value ?
          false
        :
          true
        )
      );
    }

  };


//! @}
