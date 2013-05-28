/*! 
  \file Main.h Main Doxygen documentation page
*/

//! Main namespace for So-Bogus
namespace bogus {

/*! 
\mainpage 

\section intro Introduction

So-Bogus regroups a set of loosely connected components, architectured as follow:

- \ref core Header-only library, 
  - \ref block A block-sparse matrix library
  - \ref block_solvers Solvers ( GaussSeidel, ConjugateGradient ) using those matrices
  - \ref soc Tools for solving Second Order Cone complementarity problems with the above solvers
- \ref interfaces Convenient, compiled wrappers for the most popular uses of the Core libraries, such as solving Coulomb friction problems. 

\section core Core

\ref core is a set of heavily templated, header-only libraries.
\ref block is self-contained and independant from the rest of the libraries,
\ref block_solvers makes use of \ref block, and \ref soc is probably rather useless without the two previous ones. 

\subsection core_conventions Header Naming Conventions
For each module of the \c Core library, several header files are available, follwing this naming pattern:
- \c Module.fwd.hpp Forward declarations of the module's public classes
- \c Module.hpp Public classes definitions ; does not include templated method mimplmentation
- \c Module.impl.hpp Full implementation

For the \ref block module, the file Block.io.hpp includes additional definitions for IO and serialization-related functions. 

\subsection core_configuration Configuration

\subsubsection block_solvers BlockSolvers

\subsubsection soc SecondOrder

\section interfaces Interfaces

The Interfaces part of bogus is not considered stable yet, but can be used 
as examples of how to use the \ref core libraries.



*/


} //namespace bogus


