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


//! \addtogroup running_stat
//! @{



//! Class for keeping statistics of a continuously sampled process / signal.
//! Useful if the storage of individual samples is not necessary or desired.
//! Also useful if the number of samples is not known beforehand or exceeds 
//! available memory.
template<typename eT>
class running_stat
  {
  public:
  
  typedef typename get_pod_type<eT>::pod_type T;
  
  
  inline      running_stat();
  
  inline void operator() (const T sample);
  inline void operator() (const std::complex<T>& sample);

  inline void reset();
  
  inline eT   mean() const;
  
  inline  T   var   (const u32 norm_type = 0) const;
  inline  T   stddev(const u32 norm_type = 0) const;
  
  inline eT   min()  const;
  inline eT   max()  const;
  
  //
  //
  
  private:
  
  arma_aligned  T  N;
  
  arma_aligned eT  acc1;
  arma_aligned  T  acc2;
  
  arma_aligned eT  min_val;
  arma_aligned eT  max_val;
  
  arma_aligned  T  min_val_norm;
  arma_aligned  T  max_val_norm;


  friend class running_stat_aux;
  };



class running_stat_aux
  {
  public:
  
  template<typename eT>
  inline static void update_stats(running_stat<eT>&               x,  const eT               sample);

  template<typename T>
  inline static void update_stats(running_stat< std::complex<T> >& x, const T                sample);

  template<typename T>
  inline static void update_stats(running_stat< std::complex<T> >& x, const std::complex<T>& sample);

  //

  template<typename eT>
  inline static eT var(const running_stat<eT>&                x, const u32 norm_type = 0);

  template<typename T>
  inline static  T var(const running_stat< std::complex<T> >& x, const u32 norm_type = 0);

  };



//! @}
