=== Contents ===

1. General Caveats
   1.1. Windows Specific Caveats
2. Installation
   2.1. Installation via CMake 
   2.2. Manual Installation
3. Example Programs
4. Reference Manual / User Documentation
5. Technical Documentation
6. Using Armadillo in Conjunction with IT++
7. License
8. Bug Reports
9. Credits



=== 1. General Caveats ===

While this library has gone through testing, not all possible
cases have been covered yet. As such, its functionality may not 
be 100% correct.  If you find a bug, either in the library or 
the documentation, we are very interested in hearing about it.



=== 1.1. Windows Specific Caveats ===

The development and testing has so far been mainly done on UNIX-like 
platforms (Linux and MacOS), however there should be little or no 
platform specific code.  While rudimentary tests were done on a 
Windows machine, the developers are interested in hearing how well 
Armadillo works in more thorough tests.

The 'MS Visual C++ 2008 Express Edition' compiler is known
not to work.  This is due to its incomplete implementation 
of the C++ standard.

Alternative free compilers include:
  - Intel's C++ compiler
    http://software.intel.com/en-us/intel-compilers/

  - GCC (part MingGW)
    http://www.mingw.org/

  - GCC (part of CygWin)
    http://www.cygwin.com/



=== 2. Installation ===

If you have installed Armadillo using an RPM or DEB package,
you don't need to do anything else.  Otherwise read on.

There are two ways of installing Armadillo:
  (a) via CMake
  (b) manual installation

The functionality of Armadillo will be affected by what
libraries are present on your system.  Before installing
Armadillo, it's recommended that the following libraries 
are present: LAPACK, BLAS, ATLAS and Boost.

LAPACK and BLAS are the most important. If you have ATLAS and Boost,
it's also necessary to have the corresponding header files installed.
See the "Example Programs" section for more info.



=== 2.1. Installation via CMake ===

"cmake" (version 2.6 or later) needs to be present on your system.
It's available as a pre-built package in major Linux and UNIX 
distributions, though the package may need to be explicitly installed.
If you cannot find a pre-built package for your system,
CMake can be downloaded from http://www.cmake.org

To install Armadillo, open a shell with administrator privileges
(e.g. root), change into the directory that was created by unpacking
the armadillo archive, and type the following commands:

  cmake .
  make 
  make install

CMake will figure out what other libraries are currently installed
and will modify Armadillo's configuration correspondingly.
(i.e. before installing Armadillo, install BLAS and LAPACK if you can)

If you don't have administrator privileges, change
  make install
to
  make install DESTDIR=another_location

where "another_location" is a directory where you have
write access.



=== 2.2. Manual installation ===

You will need to modify "include/armadillo_bits/config.hpp" 
to indicate which libraries are currently available on your system
and then copy the entire "include" directory to a convenient location.
Note that to use the ATLAS and Boost libraries, their header files
also need to be present.



=== 3. Example Programs ===

The "examples" directory contains several quick example programs 
that use the Armadillo library. Please see "examples/Makefile", 
which may may need to be configured for your system.

If Armadillo header files were installed in a non-standard location,
you will need to modify "examples/Makefile" to tell tell the compiler
where they are.

If Armadillo was installed manually, you will also need to explicitly 
link your programs with the libraries that were specified in 
"include/armadillo_bits/config.hpp".  For example, if you specified
that LAPACK and BLAS are available, under Linux you would use 
"-llapack -lblas" instead of "-larmadillo". Under MacOS, you would 
use "-framework Accelerate" instead of "-larmadillo".

"example1.cpp" doesn't need any external libraries.
"example2.cpp" requires the LAPACK library. You may get errors
at compile or run time if LAPACK (or its emulation via ATLAS) 
is not installed.



=== 4. Reference Manual / User Documentation ===

A quick reference manual is available at http://arma.sourceforge.net
or in the "docs_user" directory.  Use a web browser to open the 
"docs_user/index.html" file.  The documentation explains the classes
and functions, with snippets of example code. 



=== 5. Technical Documentation ===

The technical documentation (produced with the aid of Doxygen) is
available in the "docs_tech" directory.  Use a web browser to open
the "docs_tech/index.html" file.

The technical documentation helps in understanding the internals 
of Armadillo.



=== 6. Using Armadillo in Conjunction with IT++ ===

If you wish to use the IT++ library in conjunction with Armadillo,
use #include "armadillo_itpp" instead of #include "armadillo"
in your code.  See also the "examples/example_itpp.cpp" file.



=== 7. License ===

Please see the "LICENSE.txt" file.



=== 8. Bug Reports ===

If you find a bug, either in the library or the documentation,
we are very interested in hearing about it. Please send a report to:
Conrad Sanderson <conradsand at ieee dot org>



=== 9. Credits ===

Main developers:
- Conrad Sanderson, http://www.itee.uq.edu.au/~conrad/
- Ian Cullinan

Contributors:
- Justin Bedo
- Charles Gretton
- Edmund Highcock
- Kshitij Kulshreshtha 
- Oka Kurniawan
- Artem Novikov
- Martin Orlob
- Adam Piątyszek
- Ola Rinta-Koski
- Laurianne Sitbon
- Yong Kang Wong


