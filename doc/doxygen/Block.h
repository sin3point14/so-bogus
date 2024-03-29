/*!
  \file Block.h
  \brief High level documentation for the Core/Block module
*/

namespace bogus {

/*!

\page block Block
\tableofcontents

\note This module is released under the terms of the <a href="http://mozilla.org/MPL/2.0/">Mozilla Public License version 2.0</a>

\section block_introduction Introduction

The goal of this module is to make arithmetic expressions involving sparse block matrices trivial to express,
without sacrificing performance. For instance, we will be able to compute the regularized Delassus operator,
\f$ W := H M^{-1} H' + \gamma Id \f$ simply by writing
\code
W = H * ( MInv * H.transpose() ) + gamma * Id  ;
\endcode
where \c MInv is a block diagonal matrix containing LDLT factorizations
of the diagonal blocks of \c M. The evaluation will be somewhat lazy -- which
means that it will try to use as little as possible temporary storage -- and
will take advantage of parallel architectures
-- provided OpenMP support is enabled, see \ref core_configuration.

\section block_basics Basics

To use this library,
\code
#include <bogus/Core/Block.impl.hpp>

// If you need serialization capabilities
#include <bogus/Core/Block.io.hpp>
\endcode


The main user-fronting class of this library is SparseBlockMatrix, even if most of its methods are actually implemented by its parent, SparseBlockMatrixBase. FlatSparseBlockMatrix provides a more efficient implemtation with a similar interface when the blocks are dynamically-sized matrices.
Alternatively, the MappedSparseBlockMatrix class allows an externally constructed sparse block matrix to be used inside bogus expressions.

SparseBlockMatrix, FlatSparseBlockMatrix, and MappedSparseBlockMatrix are templated with a block type, and an optional set of \ref flags.
The library has been written to mainly accept Eigen dense and sparse matrices as block types, though little effort is required to make it compatible with other types.
Additionally, scalar types ( such as \c double, \c float or \c int ) are supported out of the box, as well as LU and LDLT factorizations of Eigen matrices.

The possible values for the \ref flags are a combination of:
 - \ref flags::COL_MAJOR
 - \ref flags::SYMMETRIC
 - \ref flags::UNCOMPRESSED ( except for MappedSparseBlockMatrix which maps to a compressed index format )


\section block_create Creating a SparseBlockMatrix

Three steps are necessary to create a block sparse matrix:
 - Define the block dimensions using SparseBlockMatrixBase::setRows() and SparseBlockMatrixBase::setCols()
 - Optionally pre-allocate the necessary memory using SparseBlockMatrixBase::reserve()
 - Insert the non-zeros block using SparseBlockMatrixBase::insertBack(), SparseBlockMatrixBase::insert() or one of their variations
 - <b> Call SparseBlockMatrixBase::finalize() </b>. Really. This would have been a legitimate use for the late HTML blink tag.

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
Furthermore, if you don't want the resize() method to be called at insert time, you can just use SparseBlockMatrixBase::insertBack().
The following code would work just as well, despite the increased verbosity:
\code
bogus::SparseBlockMatrix< Eigen::MatrixXd > sbm ;

// Insert the diagonal blocks
sbm.insertBack(0,0) ;
sbm.insertBack(1,1) ;

sbm.finalize() ; // Compulsory !

sbm.block(0,0) = Eigen::MatrixXd::Ones( 3, 2 ) ;
sbm.diagonal(1) = Eigen::MatrixXd::Ones( 2, 2 ) ;

\endcode

The SparseBlockMatrixBase::insertBack() insertion method requires the blocks to be inserted <b>in order</b>.
That is, for a row-major block matrix, filling the rows in increasing order, and filling each row with stricly increasing column indices.

Alternatively, <b> for matrices with an uncompressed index</b>, the SparseBlockMatrixBase::insert() method may be used.
In this case, there are no restrictions on the order in which blocks are inserted.
However, this is at the cost of poorer performance. Please refere to those methods documentation for more details.

\section block_flat ``Flat'' SparseBlockMatrix
Internally, the SparseBlockMatrix<BlockT> class uses a \c std::vector<BlockT>
to store its data. This work great for fixed-size blocks, as the storage is then contiguous.
However, dynamically-sized blocks will be individually allocated on the heap, which may 
degrade the caching performance when performing operations such as matrix-vector multiplications.

The class FlatSparseBlockMatrix<BlockT> remedies to this problem. Its public interface is similar to 
that of SparseBlockMatrix, but internally it uses a contiguous storage even for dynamically sized matrices, which may improve performance in this case. 
Note that using FlatSparseBlockMatrix to store fixed-size blocks is wasteful, as it involves
superfluous arithmetic operations keeping track of the individual blocks dimensions.
SparseBlockMatrix and FlatSparseBlockMatrix should be interchangeable in most cases 
(e.g., the result of an addition of SparseBlockMatrix can be affected to a FlatSparseBlockMatrix, and vice versa). However, FlatSparseBlockMatrix is much more experimental at this point. 
Moreover, while SparseBlockMatrix can accomodate arbitrary block types (such as other sparse matrices, or references to linear solvers), the blocks inside a FlatSparseBlockMatrix must be simple dense matrices, such as Eigen::MatrixXd.


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

\section block_access Reading the contents of a SparseBlockMatrix

Iterators over the inner dimension of a SparseBlockMatrix can be constructed
using the SparseBlockMatrixBase::innerIterator() method. For instance,
the contents of a row-major matrix can be read as follow

\code
typedef bogus::SparseBlockMatrix< Eigen::MatrixXd > SBM ;
SBM sbm ;

//...

for( int row = 0 ; row < sbm.rowsOfBlocks() ; ++ row ) {
	for( SBM::InnerIterator it ( sbm.innerIterator( row ) ) ; it ; ++it ) {
		std::cout << "Block at (" << row << ", " << it.inner() << ") is \n "
				  << sbm.block( it.ptr() ) << "\n" ;
	}
}

\endcode

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

Any SparseBlockMatrix can be assigned to another one (or to a FlatSparseBlockMatrix), as long as their block types are compatible.

\code
bogus::SparseBlockMatrix< Eigen::Matrix3d, bogus::UNCOMPRESSED  > sbm ;
// Fill sbm
// ...

// Assignment
bogus::SparseBlockMatrix< Eigen::MatrixXd, bogus::COL_MAJOR >
  sbm2 = sbm ;

\endcode

\subsection block_scaling Coefficient-wise scaling

The multiplication of each block of a SparseBlockMatrix with a scalar can be conveniently done using the SparseBlockMatrixBase::scale() method or the \c '*=' and \c '/=' operators.


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

Multiplication with a Eigen dense vector ( or matrix ) can ben done using simply the standard \c '*' operator,
or using the BlockMatrixBase::multiply() method for more flexibility ( BLAS ??mv-like ) and control over
temporary memory allocation.

Some examples
\code

Eigen::VectorXd res  = sbm * Eigen::VectorXd::Ones( n ) ;
Eigen::VectorXd res2 = sbm.transpose() * Eigen::VectorXd::Map( data, n )

sbm.multiply< true >( res2, res, -1, 1 ) ; // res -= sbm.transpose() * res2

\endcode

\note When using the \c * operator, the result of the operation will
be evaluated lazily, that is not until it is assigned to another vector or
evaluated as part of a larger expression. However, for aliasing safety reasons and
consistency with the Eigen library, the matrix-vector product may sometimes
be unecessarily be evaluated inside a temporary. You can change this behavior
using the noalias() operator of Eigen lvalues:
\code
Eigen::VectorXd rhs, res ;
// Equivalent to Eigen::VectorXd tmp = sbm*rhs ; res += sbm ;
res += sbm * rhs ;

// Equivalent to sbm.multiply< false >( rhs, res, 1, 1 ) ;
res.noalias() += sbm * rhs ;
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

\subsection block_expressions Expressions
bogus uses lazy evaluation, and as such does not evaluate operations between matrices immediately, but stores them inside expression templates.
While these expressions are transparently resolved at assignment time, and therefore should not generally matter to the end-user, they may be useful for the definition of matrix-free iterative solvers (see e.g. \ref block_solvers_is).

bogus defines two unary operations, \ref Transpose and \ref Scaling, and two binary operations, \ref Product and \ref Addition. For instance,
\code 
Addition< Scaling<AType>, Product< BType, Transpose<CType> > > expr = 3*A + B*C.transpose() ;
\endcode

In certain situations, one may not know at compile time the number of
operations that define an expression. For such scenarios, bogus defines the \ref NarySum expression. For instance, if we want to use to use the expression
\f$ H H^T + \sum_i{ a_i J_i M_i J_i^T } \f$ as a system matrix, we can do
\code
// Common type for each of the sum's operands
typedef bogus::Product< JType,
                  bogus::Product< MType, bogus::Transpose< JType > > >
                  JMJtProd ;

// Construct n-ary sum expression
// The (common) number of rows and columns of the operands has to be provided beforehand, in order to allow empty sum expressions
bogus::NarySum< JMJtProd > sum( nRows, nCols ) ;
for( unsigned i = 0 ; i < Jmatrices.size() ; ++i ) {
    const JType &J = Jmatrices[i] ;
    sum += a[i] * ( J * ( Mmatrices[i] * J.transpose() ) );
}
  
// Construct global expression
typedef bogus::Product< HType, bogus::Transpose< HType > > HHtProd ;
typedef bogus::Addition< HHtProd, bogus::NarySum< JMJtProd > > Expr ;
  
const Expr W = H * H.transpose() + sum ;

// Use the expression inside a linear solver
bogus::Krylov<Expr>(W).solve(b,x) ;

\endcode

\subsection block_limitations Limitations

As a rule of thumb, these limitations can be circumvented by explicitely assigning the result of each operation to a temporary object.

\subsubsection block_limit_aliasing Aliasing
For performance reasons, operations on SparseBlockMatrix should not be assumed to be aliasing-safe.
This is especially true for matrix-vector multiplication using directly BlockMatrixBase::multiply()
( rhs and res should not alias ), and matrix-matrix addition
( the matrix that is being assigned to should not appear anywhere but as the left-most operand ).

\subsubsection block_limit_type_deduction Type deduction
In some cases bogus will not be able to deduce the correct return type for an operation.
This can happen for matrix-matrix products that have to be evaluated as a part of a larger arithmetic expressions,
and which involve "unusual" block types. If such an error occurs, just assign the offending product to a temporary SparseBlockMatrix.

\subsubsection block_limit_performance Performance
bogus will not necessarily choose the most optimized type for the evaluation of temporary expressions. For example,
it might choose a row-major matrix when a column-major one would be more approriate, or fail to notice that
\c H \c * \c H.transpose() should use symmetric storage.

Explicit parenthesisation will also help performance. Otherwise, bogus may 
perform matrix-matrix products in an inefficient order, or allocate
unnecessary temporary matrices.

\section block_objects Other useful objects


\subsection block_iterable_objects Matrix-like objects
The Zero class is a useful placeholder for matrix arguments that are not always useful.
For instance, in the GaussSeidel solver, problems with and without linear constraints use
the same code; in the latter case, the constraint matrix is represented with the Zero class.

CompoundBlockMatrix is a utility class representing the concatenation of two objects, which
may have different block types.

*/



} //namespace bogus

