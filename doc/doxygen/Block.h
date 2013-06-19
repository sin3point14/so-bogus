/*! 
  \file Block.h 
  \brief High level documentation for the Core/Block module
*/

namespace bogus {

/*! 

\page block Block
\tableofcontents

\section block_basics Basics

To use the library, 
\code
#include <bogus/Core/Block.impl.hpp>

// If you need serialization capabilities
#include <bogus/Core/Block.io.hpp>
\endcode


The main user-fronting class of this library is SparseBlockMatrix, even if most its method are actually implemented by SparseBlockMatrixBase.

SparseBlockMatrix are templated with a block type, and an optional set of \ref flags. The library has been written to mainly accept Eigen dense and sparse matrices as block types, though little effort is required to make it compatible with other types. Additionally, scalar types ( such as \c double, \c float or \c int ) are supported out of the box, as well as LU and LDLT factorizations of Eigen matrices.

The possible values for the \ref flags are a combination of:
 - \ref flags::COMPRESSED 
 - \ref flags::COL_MAJOR
 - \ref flags::SYMMETRIC


\section block_create Creating a SparseBlockMatrix

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

\section block_operations Using a SparseBlockMatrix

Currently, the following operations are supported:
 - \ref block_assign
 - \ref block_scaling
 - \ref block_transpose
 - \ref block_add
 - \ref block_mv
 - \ref block_mm

\subsection block_assign Assignment

Any SparseBlockmatrix can be assigned to another one, as long as their block types are compatible.

\code
bogus::SparseBlockMatrix< Eigen::Matrix3d > sbm ;
// Fill sbm 
// ...

// Assignment
bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::COMPRESSED | bogus::flags::COL_MAJOR > 
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

bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::flags::SYMMETRIC > sbm2 ;

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
bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::flags::SYMMETRIC | bogus::flags::COMPRESSED > W = H * H.transpose() ;
\endcode

\warning At the moment, this operation suffer from several limitations:
 - Aliasing is not handled yet ( A = A*B probably won't do what you expect )
 - Multiplications cannot be chained yet; temporary matrices much be manually created and assigned each time 
  ( A = B*C*D won't compile, but T = B*C ; A = T*D will ) 


*/


} //namespace bogus

