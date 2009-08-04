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


template<typename T1>
inline
podarray<T1>::~podarray()
  {
  arma_extra_debug_sigprint_this(this);
  
  if(n_elem > sizeof(mem_local)/sizeof(T1) )
    {
    delete [] mem;
    }

  if(arma_config::debug == true)
    {
    access::rw(mem) = 0;
    }
  }



template<typename T1>
inline
podarray<T1>::podarray()
  : n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint_this(this);
  }
  
  

template<typename T1>
inline
podarray<T1>::podarray(const podarray& x)
  : n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint();
  
  this->operator=(x);
  }
  
  
  
template<typename T1>
inline
const podarray<T1>&
podarray<T1>::operator=(const podarray& x)
  {
  arma_extra_debug_sigprint();
  
  if(this != &x)
    {
    init(x.n_elem);
    
    for(u32 i=0; i<n_elem; ++i)
      {
      access::rw(mem[i]) = x.mem[i];
      }
    }
        
  return *this;
  }
  


template<typename T1>
arma_inline
podarray<T1>::podarray(const u32 new_n_elem)
  : n_elem(0)
  , mem(0)
  {
  arma_extra_debug_sigprint_this(this);
  
  init(new_n_elem);
  }


template<typename T1>
arma_inline
T1
podarray<T1>::operator[] (const u32 i) const
  {
  return mem[i];
  }



template<typename T1>
arma_inline
T1&
podarray<T1>::operator[] (const u32 i)
  {
  return access::rw(mem[i]);
  }
  
  
  
template<typename T1>
arma_inline
T1
podarray<T1>::operator() (const u32 i) const
  {
  arma_debug_check( (i >= n_elem), "podarray::operator(): index out of bounds");
  return mem[i];
  }



template<typename T1>
arma_inline
T1&
podarray<T1>::operator() (const u32 i)
  {
  arma_debug_check( (i >= n_elem), "podarray::operator(): index out of bounds");
  return access::rw(mem[i]);
  }



template<typename T1>
inline
void
podarray<T1>::set_size(const u32 new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  init(new_n_elem);
  }



template<typename T1>
inline
void
podarray<T1>::fill(const T1 val)
  {
  arma_extra_debug_sigprint();
  
  for(u32 i=0; i<n_elem; ++i)
    {
    access::rw(mem[i]) = val;
    }
  }



template<typename T1>
inline
void
podarray<T1>::zeros()
  {
  arma_extra_debug_sigprint();
  
  fill(0);
  }



template<typename T1>
inline
void
podarray<T1>::zeros(const u32 new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  init(new_n_elem);
  fill(0);
  }



template<typename T1>
arma_inline
T1*
podarray<T1>::memptr()
  {
  return const_cast<T1*>(mem);
  }
  
  

template<typename T1>
arma_inline
const T1*
podarray<T1>::memptr() const
  {
  return mem;
  }



template<typename T1>
inline
void
podarray<T1>::init(const u32 new_n_elem)
  {
  arma_extra_debug_sigprint();
  
  if(n_elem == new_n_elem)
    {
    return;
    }
    
  if(n_elem > sizeof(mem_local)/sizeof(T1) )
    {
    delete [] mem;
    }
  
  if(new_n_elem <= sizeof(mem_local)/sizeof(T1) )
    {
    access::rw(mem) = mem_local;
    }
  else
    {
    access::rw(mem) = new T1[new_n_elem];
    }
  
  access::rw(n_elem) = new_n_elem;
  
  
  }

//! @}
