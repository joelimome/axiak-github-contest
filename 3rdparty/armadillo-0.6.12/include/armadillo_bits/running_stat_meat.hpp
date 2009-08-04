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



template<typename eT>
running_stat<eT>::running_stat()
  : N           (typename running_stat<eT>::T(0))
  , acc1        (                          eT(0))
  , acc2        (typename running_stat<eT>::T(0))
  , min_val     (                          eT(0))
  , max_val     (                          eT(0))
  , min_val_norm(typename running_stat<eT>::T(0))
  , max_val_norm(typename running_stat<eT>::T(0))
  {
  arma_extra_debug_sigprint_this(this);
  }



//! update statistics to reflect new sample
template<typename eT>
inline
void
running_stat<eT>::operator() (const typename running_stat<eT>::T sample)
  {
  arma_extra_debug_sigprint();

  running_stat_aux::update_stats(*this, sample);
  }



//! update statistics to reflect new sample (version for complex numbers)
template<typename eT>
inline
void
running_stat<eT>::operator() (const std::complex< typename running_stat<eT>::T >& sample)
  {
  arma_extra_debug_sigprint();

  isnt_same_type<eT, std::complex< typename running_stat<eT>::T > >::check();

  running_stat_aux::update_stats(*this, sample);
  }



//! set all statistics to zero
template<typename eT>
inline
void
running_stat<eT>::reset()
  {
  arma_extra_debug_sigprint();

  typedef typename running_stat<eT>::T T;
  
  N            =  T(0);
  
  acc1         = eT(0);
  acc2         =  T(0);
  
  min_val      = eT(0);
  max_val      = eT(0);
  
  min_val_norm =  T(0);
  max_val_norm =  T(0);
  }



//! mean or average value
template<typename eT>
inline
eT
running_stat<eT>::mean() const
  {
  arma_extra_debug_sigprint();

  typedef typename running_stat<eT>::T T;
  
  if(N > T(0))
    {
    return acc1 / N;
    }
  else
    {
    return eT(0);
    }
  }



//! variance
template<typename eT>
inline
typename running_stat<eT>::T
running_stat<eT>::var(const u32 norm_type) const
  {
  arma_extra_debug_sigprint();

  return running_stat_aux::var(*this, norm_type);
  }



//! standard deviation
template<typename eT>
inline
typename running_stat<eT>::T
running_stat<eT>::stddev(const u32 norm_type) const
  {
  arma_extra_debug_sigprint();

  return std::sqrt( running_stat_aux::var(*this, norm_type) );
  }



//! minimum value
template<typename eT>
inline
eT
running_stat<eT>::min() const
  {
  arma_extra_debug_sigprint();

  return min_val;
  }



//! maximum value
template<typename eT>
inline
eT
running_stat<eT>::max() const
  {
  arma_extra_debug_sigprint();

  return max_val;
  }



//! update statistics to reflect new sample
template<typename eT>
inline
void
running_stat_aux::update_stats(running_stat<eT>& x, const eT sample)
  {
  arma_extra_debug_sigprint();

  if(x.N > eT(0))
    {
    if(sample < x.min_val)
      {
      x.min_val = sample;
      }
    
    if(sample > x.max_val)
      {
      x.max_val = sample;
      }
    }
  else
    {
    x.min_val = sample;
    x.max_val = sample;
    }
  
  x.N++;
  
  x.acc1 += sample;
  x.acc2 += sample*sample;
  }



//! update statistics to reflect new sample (version for complex numbers)
template<typename T>
inline
void
running_stat_aux::update_stats(running_stat< std::complex<T> >& x, const T sample)
  {
  arma_extra_debug_sigprint();

  running_stat_aux::update_stats(x, std::complex<T>(sample));
  }



//! alter statistics to reflect new sample (version for complex numbers)
template<typename T>
inline
void
running_stat_aux::update_stats(running_stat< std::complex<T> >& x, const std::complex<T>& sample)
  {
  arma_extra_debug_sigprint();

  const T sample_norm = std::norm(sample);
  
  if(x.N > T(0))
    {
    if(sample_norm < x.min_val_norm)
      {
      x.min_val_norm = sample_norm;
      x.min_val      = sample;
      }
    
    if(sample_norm > x.max_val_norm)
      {
      x.max_val_norm = sample_norm;
      x.max_val      = sample;
      }
    }
  else
    {
    x.min_val      = sample;
    x.max_val      = sample;
    x.min_val_norm = sample_norm;
    x.max_val_norm = sample_norm;
    }
  
  x.N++;
  
  x.acc1 += sample;
  x.acc2 += sample_norm;
  }



//! variance
template<typename eT>
inline
eT
running_stat_aux::var(const running_stat<eT>& x, const u32 norm_type)
  {
  arma_extra_debug_sigprint();

  if(x.N > eT(1))
    {
    const eT norm_val = (norm_type == 0) ? (x.N-eT(1)) : x.N;
    const eT var_val  = (x.acc2 - x.acc1*x.acc1/x.N) / norm_val;
    return var_val;
    }
  else
    {
    return eT(0);
    }
  }



//! variance (version for complex numbers)
template<typename T>
inline
T
running_stat_aux::var(const running_stat< std::complex<T> >& x, const u32 norm_type)
  {
  arma_extra_debug_sigprint();

  if(x.N > T(1))
    {
    const T norm_val = (norm_type == 0) ? (x.N-T(1)) : x.N;
    const T var_val  = (x.acc2 - std::norm(x.acc1)/x.N) / norm_val;
    return var_val;
    }
  else
    {
    return T(0);
    }
  }



//! @}
