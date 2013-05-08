
#if EIGEN_VERSION_AT_LEAST(3,1,0)
#include <Eigen/SparseCore>
#else
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <Eigen/Sparse>
#endif
