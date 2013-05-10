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
\ref block_solvers makes use of \ref blocks, and \ref soc is probably rather useless without the two previous ones. 

\subsection block Block

The main user-fronting class of this library is SparseBlockMatrix, even if most oth its method are actually implemented by SparseBlockMatrixBase.

SparseBlockMatrix are templated with a block type, and an optional set of \ref flags. The library has been written to mainly accept Eigen dense and sparse matrices as block types, though little effort is required to make it compatible with other types. Aditionally, scalar types ( such as \c double, \c float or \c int ) are supported, as well as LU and LDLT factorizations.


\code

\endcode

\subsection block_solvers BlockSolvers

\subsection soc SecondOrder

\section interfaces Interfaces

The Interfaces part of bogus is not considered stable yet, but can be used 
as examples of how to use the \ref core libraries.



*/


} //namespace bogus


