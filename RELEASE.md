## Version 1.3.0
July 21, 2016

### Notable new features
 - ProductGaussSeidel avoids explicitely assembling the Delassus operator
 - GaussSeidel can now handle systems with linear equality constraints
 - Experimental ADMM and DualAMA solvers
 - ProjectedGradient may be used to solve non-optimization problems (risky)

### Bugfixes
 - Fix bug when using infinity norm without constraints

## Version 1.2.0
December 31, 2015

### Notable new features
 - Compatibility wih upcoming Eigen 3.3 branch
 - Matrix-free Projected-Gradient and iterative linear solvers 
 - New Projected-Gradient variants
 - Extracted interfaces for Cadoux algorithm 
 - Added LCP local solver 

### Bugfixes
 - Fixed a compile error with fixed-size-blocks matrices assignment
 - Do not assume existence of diagonal blocks in BlockSolvers
 - Fixed nested evaluation of matrix-matrix products with symmetric lhs
 - Improved compatibility between matrices with block of types `Scalar` and `Eigen::Matrix<Scalar,1,1>` 

## Version 1.1.0
August 8, 2015

### Notable new features
 - Added several variants of ProjectedGradient algorithm
 - Compatibility with `fclib` 0.2

### Bugfixes
 - Fixed compilation with `clang -std=c++11`
 - Fixed conversion from `bogus::SparseBlockMatrix` to `Eigen::SparseMatrix`

## Version 1.0.0
December 13, 2013
