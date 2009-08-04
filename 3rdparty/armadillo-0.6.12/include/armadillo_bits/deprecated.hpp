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


//! \addtogroup deprecated
//! @{




#define basic_mat    Mat
#define basic_colvec Col
#define basic_rowvec Row
#define basic_math   Math


namespace deprecated
  {
  // pre-processing tricks to let doxygen document deprecated classes.
  // the defines above were not documented via \def as we've explicitly
  // set Doxygen's ENABLE_PREPROCESSING to NO
  
  #if defined(ARMA_TEMP_DEFINE)
    #undef ARMA_TEMP_DEFINE
  #endif
  
  #if defined(ARMA_TEMP_DEFINE)
    
    #undef basic_mat
    #undef basic_colvec
    #undef basic_rowvec
    
    //! Class name 'basic_mat' is deprecated.
    //! 'basic_mat' has been defined as 'Mat' for compatibility with Armadillo <= 0.6.2
    template<typename eT> class basic_mat {};
    
    //! Class name 'basic_colvec' is deprecated.
    //! 'basic_colvec' has been defined as 'Col' for compatibility with Armadillo <= 0.6.2
    template<typename eT> class basic_colvec {};
    
    //! Class name 'basic_rowvec' is deprecated.
    //! 'basic_rowvec' has been defined as 'Row' for compatibility with Armadillo <= 0.6.2
    template<typename eT> class basic_rowvec {};
  
    //! Class name 'basic_math' is deprecated.
    //! 'basic_math' has been defined as 'Math' for compatibility with Armadillo <= 0.6.2
    template<typename eT> class basic_math {};
  
  
  #endif
  }


//! @}
