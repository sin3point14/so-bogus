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

The main user-fronting class of this library is SparseBlockMatrix, even if most its method are actually implemented by SparseBlockMatrixBase.

SparseBlockMatrix are templated with a block type, and an optional set of \ref flags. The library has been written to mainly accept Eigen dense and sparse matrices as block types, though little effort is required to make it compatible with other types. Additionally, scalar types ( such as \c double, \c float or \c int ) are supported out of the box, as well as LU and LDLT factorizations of Eigen matrices.

The possible values for the \ref flags are a combination of:
 - \ref flags::COMPRESSED 
 - \ref flags::COL_MAJOR
 - \ref flags::SYMMETRIC

\subsubsection block_create Creating a SparseBlockMatrix

Three steps are necessary to create a block sparse matrix:
 - Define the block dimensions use SparseBlockMatrixBase::setRows() and SparseBlockMatrixBase::setCols()
 - Insert the non-zeros block using SparseBlockMatrixBase::insertBack() or SparseBlockMatrixBase::insertBackAndResize()
 - Call SparseBlockMatrixBase::finalize() 

Example: Creating a block-diagonal matrix  
\code

//  Creates the following matrix:
//
//  1 1 0 0  
//  1 1 0 0
//  1 1 0 0
//  0 0 1 1
//  0 0 1 1
//

int rowsPerBlock[2] = { 3, 2 } ;

bogus::SparseBlockMatrix< Eigen::MatrixXd > sbm ;
sbm.setRows( 2, rowsPerBlock ) ; // 2 block-rows which each have a number of rows given by rowsPerBlock
sbm.setCols( 2, 2 ) ; // 2 block-columns which all have 2 columns

// Insert the diagonal blocks
sbm.insertBackAndResize( 0, 0 ).setOnes() ;
sbm.insertBackAndResize( 1, 1 ).setOnes() ;

sbm.finalize() ; // Compulsory !

\endcode

Note that you don't have to specify the block contents before the call to SparseBlockMatrixBase::finalize(). 
Furthermore, if you don't want the \c resize() method of the blocks to be called at insert time, you can just use SparseBlockMatrixBase::insertBack(). The following code would work just as well, despite the increased verbosity:
\code 
bogus::SparseBlockMatrix< Eigen::MatrixXd > sbm ;

// Insert the diagonal blocks
sbm.insertBack(0,0) ;
sbm.insertBack(1,1) ;

sbm.finalize() ; // Compulsory !

sbm.block(0,0) = Eigen::MatrixXd::Ones( 3, 2 ) ;
sbm.diagonal(1) = Eigen::MatrixXd::Ones( 2, 2 ) ;

\endcode


The insertion method is named insertBack() to warn you that in the general case, the blocks have to be inserted in order.
That is, for a row-major block matrix, filling the rows in increasing order, and filling each row with stricly increasing column indices.
Note that this constraint can be safely ignored when using an uncompressed index. ( See \ref flags and SparseBlockMatrix::insertBack() ) 

\subsubsection block_operations Using a SparseBlockMatrix

Currently, the following operations are supported:
 - Assignment
 - Transposition
 - Block Matrix/Dense Vector multiplication
 - Block Matrix/Block Matrix multiplication

\b Assignment \n

Any SparseBlockmatrix can be assigned to another one, as long as their block types are compatible.
\warning Assigning a SparseBlockMatrix with \ref flags::SYMMETRIC to a non-symmetric one will result in only the triangular part being copied, instead of the full matrix

\code
bogus::SparseBlockMatrix< Eigen::Matrix3d > sbm ;
// Fill sbm 
// ...

// Assignment
bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COMPRESSED | bogus::flags::COL_MAJOR > 
  sbm2 = sbm ;

\endcode

\b Transposition

Transposing a block matrix can be done with the SparseBlockMatrixBase::transpose() method.
Note that this method does not do any work; it can only be use as the right hand side of an affectation operation,
or inside an arithmetic operation ( see below )

\b Block \b Matrix/Dense \b Vector \b multiplication

Multiplication against a Eigen dense vector can ben done using simply the standard \c * operator,
or using the BlockMatrixBase::multiply() method for more flexibility and no temporary memory allocation.

Some examples
\code

Eigen::VectorXd res  = sbm * Eigen::VectorXd::Ones( n ) ;
Eigen::VectorXd res2 = sbm.transpose() * Eigen::VectorXd::Map( data, n )

sbm.multiply< true >( res2, res, -1, 1 ) ; // res -= sbm.transpose() * res2

\endcode


\subsection block_solvers BlockSolvers

\subsection soc SecondOrder

\section interfaces Interfaces

The Interfaces part of bogus is not considered stable yet, but can be used 
as examples of how to use the \ref core libraries.



*/


} //namespace bogus


