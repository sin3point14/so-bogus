## Compiling and Installing

### Header only library
Obviously, no compilation is required if you only plan to use the header-only part of bogus.
Using the CMake build script can still be helpful for copying all the header files to a target
install directory, as detailed in the next section.

A `bogus -> src` symlink is provided for convenience, so you can `#include <bogus/Core/Block.impl.hpp>`
on both local and `make install`ed copies of bogus.

### CMake build script

The CMake build script serves several purposes:

 - Copy and install the header-only library ( `src/Core` and `src/Extra` folders )
 - Build and install the So-bogus library ( `src/Interfaces` folder )
 - Build the test scripts ( `tests/` folder )
 - Build the standalone applications ( `apps/` folder )

To use this script, simply create an empty `build` directory, then run
 `cmake /path/to/bogus [OPTIONS]`, `make`, and optionally `make install`.

#### Dependencies

The script will attempt to find the following libraries:

 - [Eigen 3.x] (http://eigen.tuxfamily.org/) **(REQUIRED)**
   
 If not installed in a standard location, you can define the `EIGEN_ROOT` environment variable.
 - __OpenMP__ (optional)

 If the C++ compiler supports OpenMP, the So-bogus library will be built with parallelization enabled.
 This can be overriden with the `-DOPENMP=off` cmake flag.
 - [boost\_serialization](http://www.boost.org/) (optional)

  The `boost_serialization` library is required for writing and reading bogus structures from files.  This can be disabled with the `-DBOOST_SERIALIZATION=off` cmake flag.
 - __Intel® MKL__ (optional)

   Sparse block matrices that follow the Block Sparse Row (BSR) format may use MKL
   for some operations, such as matrix-vector product. This can be disabled with the `-DMKL=off` cmake flag. If the MKL installation directory is not standard, it can be specified with the `MKL_ROOT` environment variable.
 - [FCLib](http://fclib.gforge.inria.fr/) (optional)

  The `apps/FCLibLoader` program can read and solve problems from the FCLib 
  library of 3D frictional contact problems.

If the tests are enabled, the build script will also download and compile the [googletest](https://code.google.com/p/googletest/) framework. A working Internet connection is therefore required.
   
#### Options
  
In addition to the aforementioned options that disable some optional dependencies, 
the following CMake options are available. The can be defined using the  command line syntax `-DOPTION_NAME=value`.

 - __LIB__=(__on__|off)  Whether to build the So-bogus dynamic library (`libbogus.so`).
 - __TESTS__=(__on__|off)  Whether to build the somewhat unit tests.
 - __APPS__=(__on__|off)  Whether to build the standalone applications. You probably won't need them.
 - __WITH_2D__=(__on__|off) Compile support for 2D problems  in the So-bogus library
 - __WITH_3D__=(__on__|off) Compile support for 3D problems support in the So-bogus library
 - __WITH_DYNAMIC__=(on|__off__) Compile support for dynamically-sized problems in the So-bogus library
 
A few compiler flags (`gcc` and `clang`) can be set using the following options:

 - __CPP11__=(on|__off__) Compile with `-std=c++11` 
 - __STRICT__=(on|__off__) Enable more warnings and treat them as errors
 - __FAST_MATH__=(__on__|off) Enable optimizations that break IEEE 754 compliance (`-ffast-math`) 

You may also find the following standard CMake variables useful:

 - __CMAKE\_BUILD\_TYPE__=(__Release__|Debug) The Release mode enables compiler optimizations and disables runtime assertions.
 - __CMAKE\_INSTALL\_PREFIX__=/path Defines the path to which files will be copied when using the `make install` command. More precisely, the header files from the `src/` directory will be copied to `/path/include/`, and the dynamic `libbogus.so` library to `/path/lib/`.

### Testing
You can check that bogus works relatively well on your machine by running `tests/testbogus`. Any failing test report would be appreciated.

