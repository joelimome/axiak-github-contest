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


//! \addtogroup podarray
//! @{



//! A lightweight array for POD types. If the amount of memory requested is small, the stack is used.

template<typename T1>
class podarray
  {
  public:
  
  //! number of elements held
  const u32 n_elem;    

  //! pointer to memory used by the object
  arma_aligned const T1* const mem;
  
  protected:
  //! Internal memory, to avoid calling the 'new' operator for small amounts of memory.
  arma_aligned T1 mem_local[ 16 ];
  
  
  public:
  
  inline ~podarray();
  inline  podarray();
  
  inline                 podarray (const podarray& x);
  inline const podarray& operator=(const podarray& x);
  
  arma_inline explicit podarray(const u32 new_N);

  arma_inline T1& operator[] (const u32 i);
  arma_inline T1  operator[] (const u32 i) const;
  
  arma_inline T1& operator() (const u32 i);
  arma_inline T1  operator() (const u32 i) const;

  inline void set_size(const u32 new_n_elem);
  inline void fill(const T1 val);

  inline void zeros();
  inline void zeros(const u32 new_n_elem);

  arma_inline       T1* memptr();
  arma_inline const T1* memptr() const;
  
  
  protected:
  
  inline void init(const u32 new_n_elem);

  };

//! @}
