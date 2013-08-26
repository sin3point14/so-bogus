/*! 
  \file Block.h 
  \brief High level documentation for the Core/Block module
*/

namespace bogus {

/*! 

\page block Block
\tableofcontents

\section block_basics Basics

To use this library, 
\code
#include <bogus/Core/Block.impl.hpp>

// If you need serialization capabilities
#include <bogus/Core/Block.io.hpp>
\endcode


The main user-fronting class of this library is SparseBlockMatrix, even if most its method are actually implemented by SparseBlockMatrixBase. Alternatively, the MappedSparseBlockMatrix class allows an externally constructed sparse block matrix to be used inside bogus expressions.

SparseBlockMatrix and MappedSparseBlockMatrix are templated with a block type, and an optional set of \ref flags. The library has been written to mainly accept Eigen dense and sparse matrices as block types, though little effort is required to make it compatible with other types. Additionally, scalar types ( such as \c double, \c float or \c int ) are supported out of the box, as well as LU and LDLT factorizations of Eigen matrices.

The possible values for the \ref flags are a combination of:
 - \ref flags::COL_MAJOR
 - \ref flags::SYMMETRIC
 - \ref flags::UNCOMPRESSED ( except for MappedSparseBlockMatrix which maps to a compressed index format )


\section block_create Creating a SparseBlockMatrix

Three steps are necessary to create a block sparse matrix:
 - Define the block dimensions using SparseBlockMatrixBase::setRows() and SparseBlockMatrixBase::setCols()
 - Optionally pre-allocate the necessary memory using SparseBlockMatrixBase::reserve() 
 - Insert the non-zeros block using SparseBlockMatrixBase::insertBack() or SparseBlockMatrixBase::insertBackAndResize()
 - <b> Call SparseBlockMatrixBase::finalize() </b>

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


The insertion method is named insertBack() to warn the user that in the general case, the blocks have to be inserted in order.
That is, for a row-major block matrix, filling the rows in increasing order, and filling each row with stricly increasing column indices.
Note that this constraint is relaxed when using an uncompressed index. ( See \ref flags and SparseBlockMatrix::insertBack() ) 

\section block_map Creating a MappedSparseBlockMatrix

MappedSparseBlockMatrix are const references to either
 - another SparseBlockMatrixBase object that uses a compressed index
 - an external matrix in a BSR-like format 
    ( e.g. it can handle column-major matrices as well )

For instance, mapping to a SparseBlockMatrix
\code 
bogus::SparseBlockMatrix< Eigen::Matrix3d > source ;
//[...]  Fill source

bogus::MappedSparseBlockMatrix< Eigen::Matrix3d > map ;
map.mapTo( source ) ;

// can also be done using the BSR interface

map.cloneDimensions( source ) ;
map.mapTo( source.nBlocks(),
    source.data(),
    source.majorIndex().outerIndexPtr(),
    source.majorIndex().innerIndexPtr()
    );
\endcode
 
in the previous example, note that the call to SparseBlockMatrixBase::cloneDimensions() could be replaced by calls to SparseBlockMatrixBase::setRows() and SparseBlockMatrixBase::setCols().

\section block_operations Using a SparseBlockMatrix (or a MappedSparseBlockMatrix)

Currently, the following operations are supported:
 - \ref block_assign
 - \ref block_scaling
 - \ref block_transpose
 - \ref block_add
 - \ref block_mv
 - \ref block_mm

 Most of those operations can be done using the standard C++ operators, in a lazy way -- the resulting matrix will only be computed when assigned to another BlockMatrix. 
Since it is immutbale, a MappedSparseBlockMatrix can only be used in the righ-hand-side of those operations, never as a left-hand-side.
 
 \warning The operations should be able to be composed in an arbitrary way, but in pratice there are a few \ref block_limitations.

\subsection block_assign Assignment

Any SparseBlockmatrix can be assigned to another one, as long as their block types are compatible.

\code
bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::UNCOMPRESSED  > sbm ;
// Fill sbm 
// ...

// Assignment
bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::COL_MAJOR > 
  sbm2 = sbm ;

\endcode

\subsection block_scaling Coefficient-wise scaling

The multiplication of each block of a SparseBlockMatrix with a scalar can be conveniently done using the SparseBlockMatrixBase::scale() method or the \c '*=' and \c '+=' operators.  


\subsection block_transpose Transpose

Transposing a block matrix can be done with the SparseBlockMatrixBase::transpose() method.
Note that this method does not do any work; it can only be used as the right hand side of an affectation operation,
or inside an arithmetic operation ( see below )

\subsection block_add Block Matrix/Block Matrix addition

Two block matrix can be added together using the standard \c '+' operator, provided they have the same block structure
( that is, their rows and columns of blocks must have similar dimensions ). This operator does no do any work, but just returns an Addition expression which will be evaluated when it is assigned to another SparseBlockMatrix.

Alternatively, the SparseBlockMatrixBase::add method provides a ?axpy-like interface which allows on-the-fly scaling of the right-hand-side. \c '-', \c '+=' and \c '-=' are also defined.

Examples
\code

bogus::SparseBlockMatrix< Eigen::MatrixXd > sbm1 ;

// [...] fill sbm1

bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::SYMMETRIC > sbm2 ;

sbm2 = sbm1 + sbm1.transpose() ;
sbm1 -= sbm2 ;
sbm1.add< false >( sbm2, .5 ) ;

\endcode 

\subsection block_mv Block Matrix/Dense Vector multiplication

Multiplication against a Eigen dense vector can ben done using simply the standard \c '*' operator,
or using the BlockMatrixBase::multiply() method for more flexibility ( BLAS ??mv-like )and no temporary memory allocation.

Some examples
\code

Eigen::VectorXd res  = sbm * Eigen::VectorXd::Ones( n ) ;
Eigen::VectorXd res2 = sbm.transpose() * Eigen::VectorXd::Map( data, n )

sbm.multiply< true >( res2, res, -1, 1 ) ; // res -= sbm.transpose() * res2

\endcode
 
\subsection block_mm Block Matrix/Block Matrix multiplication

Two SparseBlockMatrix can me multiplied provided the dimensions of their blocks are compatible, using the standard \c '*' operator.
The return type of this operator is a \ref Product expression, which in itself does not perform any computation; the proper multiplication
will only be done when it is assigned to another SparseBlockMatrix.
\code
bogus::SparseBlockMatrix< GradBlockT > H ;
// [...] Fill H
bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::SYMMETRIC > W = H * H.transpose() ;
\endcode

\subsection block_limitations Limitations

As a rule of thumb, these limitations can be circumvented by explicitely assigning the result of each operation to a temporary object.

\subsubsection block_limit_aliasing Aliasing
For performance reasons, operations on SparseBlockMatrix should not be assumed to be aliasing-safe.
This is especially true for matrix-vector multiplication ( rhs and res should not alias ), and matrix-matrix addition
( the matrix that is being assigned to should not appear anywhere but as the left-most operand ).

\subsubsection block_limit_type_deduction Type deduction
In some cases bogus will not be able to deduce the correct return type for an operation. 
This can happen for matrix-matrix products that have to be evaluated as a part of a larger arithmetic expressions,
and which involve "unusual" block types. If such an error occurs, just assign the offending product to a temporary SparseBlockMatrix.

\subsubsection block_limit_performance Performance
bogus will not necessarily chose the most optimized type for the evaluation of temporary expressions. For example,
it might choose a row-major matrix when a column-major one would be more approriate, or fail to notice that
\c H \c * \c H.transpose() should use symmetric storage.

Explicit parenthesisation will also help performance. Otherwise, bogus may for instance perform a matrix-matrix product operation before a matrix-vector product, while the same result could be computed using only two matrix-vector operations. 

*/


} //namespace bogus

