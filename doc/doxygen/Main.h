/*! 
  \file Main.h 
  \brief Main Doxygen documentation page
*/

//! Main namespace for bogus and So-bogus
namespace bogus {

/*! 
\mainpage 

\section intro Introduction

So-bogus is a set of loosely connected components, organized as follow:

- \ref core, a MPL-licensed header-only library, which include
  - \ref block, A block-sparse matrix library
  - \ref block_solvers, Solvers ( Projected Gauss-Seidel, Krylov linear solvers ) using those matrices
- \ref extra, a GPL-licensed header-only library, which include
  - \ref soc Tools for solving Second Order Cone complementarity problems with \ref block_solvers
- \ref interfaces Convenient, compilable wrappers for the most popular uses of the header-only libraries, 
    such as solving Coulomb friction problems. GPL licensed.

\section core Core ( a.k.a. bogus )

- \ref core is a set of heavily templated, header-only libraries.
- \ref block is self-contained and independent from the other modules.
- \ref block_solvers is written to operate on matrices from \ref block.

\subsection core_conventions Header Naming Conventions
For each module of the \c Core library, several header files are available, following this naming pattern:
- \c Module.fwd.hpp Forward declarations of the module's public classes
- \c Module.hpp Public classes and operators definitions 
- \c Module.impl.hpp Full implementation

For the \ref block module, the file Block.io.hpp includes additional definitions for IO and serialization-related functions. 

\subsection core_configuration Configuration

The behavior of the header-only libraries may be customized using the following compile-time macros:

\c BOGUS_DONT_PARALLELIZE

If defined, bogus will not use multiple threads, even it is compiled with OpenMP enabled. 
<i> Not defined by default. </i>

\c BOGUS_DEFAULT_INDEX_TYPE 

Default integer index type, which will be used to refer to row or column indices. Should be signed
for compatibility with OpenMP 2.0. Should be set to \c MKL_INT if you plan to use the mkl bindings.
Can be overwritten using template parameters for e.g. MappedSparseBlockMatrix.
<i>Defaults to </i> \c int

\c BOGUS_DEFAULT_BLOCK_PTR

Integer type that will be used to refer to a given block within its sequential container.
<i>Defaults to </i> \c std::size_t

\c BOGUS_WITHOUT_EIGEN

Do not include Eigen bindings, or anything Eigen for this matter. Will break \c block_solvers and \ref soc.
<i> Not defined by default. </i>

\c BOGUS_BLOCK_WITHOUT_EIGEN_SPARSE

Do not include Eigen/Sparse bindings. 
<i> Not defined by default. </i>

\c BOGUS_WITH_MKL

Enable MKL bindings. Only supports matrix/vector multiplication for compressed SparseBlockMatrixBase with fixed-size blocks.
See \ref mkl.
<i> Not defined by default. </i>

\c BOGUS_WITH_BOOST_SERIALIZATION

Enable serialization API. 
<i> Not defined by default. </i>

\c BOGUS_SHARED_PTR_NS

If defined, bogus will use shared_ptr from this namespace instead of its own limited NaiveSharedPtr.
Values such as \c tr1 for C++03 or \c std for C++11 should work. 
<i> Not defined by default. </i>


\section extra Extra
The \ref extra header-only library collects all modules that were not suited for inclusion in \ref core,
for reasons as trivial as licensing issues.

At the time, the only module in \ref extra is \ref soc.

\ref extra follows the same header naming convention as \ref core; see \ref core_conventions.

\section interfaces Interfaces

The \ref Interfaces part of So-bogus is not considered stable yet, but provide examples of how to use the \ref core and \ref extra libraries. For instance, 
PrimalFrictionProblem is a relatively generic representation of a 2D or 3D friction
problem for which the mass matrix is a diagonal of dense blocks. It can be converted
to a DualFrictionProblem, which in turn can be solved using a GaussSeidel or ProjectedGradient block solver, or the Cadoux algorithm \cite ACML11 .

On the other hand, MecheFrictionProblem defines a more specific interface, but
provides more features, such as serialization or diagonal regularization.
Its public interface rely on as few custom classes as possible.

\note
So-bogus header-only parts are meant to be modular and composable; the \ref Interfaces library is not, and address a very specific use case.
Consequently, you should probably not use the compiled library in your program, but create your own wrapper over the \ref core and \ref extra modules. 
For instance, for maximal coordinates models,
PrimalFrictionProblem will not be optimal, as its mass matrix uses dense diagonal blocks. You will probably want to use a variation that use Eigen::SparseMatrix blocks for M and H.

Similarly, while bogus provides all the components required for solving Mixed Complementarity Problems, no such solver is exposed through the \ref Interfaces .

*/


} //namespace bogus


