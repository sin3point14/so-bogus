#include "SOCLaw.impl.hpp"

namespace bogus
{

template struct  SOCLaw< Eigen::Matrix2d,  true >  ;
template struct  SOCLaw< Eigen::Matrix3d,  true >  ;
template struct  SOCLaw< Eigen::Matrix2d, false >  ;
template struct  SOCLaw< Eigen::Matrix3d, false >  ;

}
