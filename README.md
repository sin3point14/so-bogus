###Overview

This free software consist of two complementary libraries, the first being fairly generic while the other targets more specific problems. Both are open-source and can be downloaded from the [BitBucket repository](https://bitbucket.org/gdaviet/so-bogus/downloads#tag-downloads). 

####bogus
 is a template C++ 98 header-only **sparse block matrix** library. It allows to express in a concise way algebraic operations on sparse matrices whose elements can be a variety of block types, such as [Eigen](http://eigen.tuxfamily.org) dense or sparse matrices.

####So-bogus
 is a library that uses and completes bogus to solve systems with **second order cone constraints**, such as **Coulomb friction problems**.

###Main features

#### bogus

 - Simple creation of sparse block matrices with several storage options (row or column major, symmetric, ...)
 - Compatibility with standard Block Sparse Row storage
 - Complex arithmetic expressions using standard C++ operators
 - Lazy evaluation
 - Parallelization of most sparse block matrix operations using OpenMP (if available)
 - Optional serialization using `boost::serialization`
 - Several linear solvers ( Conjugate Gradient and variations, GMRES, ...)  
 - Generic Projected Gauss-Seidel and Projected Gradient solvers for
constrained problems

#### So-bogus

 - Second Order Cones complementarity function with local solvers
 - Global Coulomb friction and SOCQP solvers

###License
**bogus** is released under the terms of  the  [Mozilla Public License v2.0](http://gdaviet.fr/code/bogus/src/master/MPL-LICENSE-2.0.txt).

**So-bogus** is released under the terms of the [GNU General Public License version 2](http://gdaviet.fr/code/bogus/src/master/GPL-LICENSE-2.txt) or, at your option, any later version.

For more information, see [LICENCE.md](http://gdaviet.fr/code/bogus/src/master/LICENSE.md).

###Documentation

**bogus** includes documentation in the doxygen format, along with usage samples.
Run `doxygen` from the `doc/doxygen` folder, or [browse it online](http://gdaviet.fr/doc/bogus/master/doxygen/).

###Requirements

Apart from a C++98 compiler, the only requirement is  [Eigen](http://eigen.tuxfamily.org).

**bogus** has been successfully tested on a variety of Linux distributions  with gcc 4.5+,
and on Mac OSX 10.8 with both gcc 4.2 and clang 3.2.
Visual Studio 2012 ( and maybe earlier ) should also be able to compile and run the test suite, though no particular attention is given to this platform.

For more information, see [INSTALL.md](http://gdaviet.fr/code/bogus/src/master/INSTALL.md).


###About
**bogus** and **So-bogus** are Copyright 2013 
[Gilles Daviet](http://gdaviet.fr) <gdaviet@gmail.com>.

**So-bogus** implements some algorithms initially developed within the [BiPop](http://bipop.inrialpes.fr) team at [Inria Rhônes-Alpes](http://inria.fr/en/centre/grenoble). Its name is derived from the main use case of the library, a **B**l**o**ck **G**a**u**ss-**S**eidel solver for **S**econd **O**rder cone problems. 




