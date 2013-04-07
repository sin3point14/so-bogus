#include "LocalSOCSolver.impl.hpp"

namespace bogus
{

template struct LocalSOCSolver< 2u, double, false > ;
template struct LocalSOCSolver< 3u, double, false > ;
template struct LocalSOCSolver< 2u, double, true  > ;
template struct LocalSOCSolver< 3u, double, true  > ;

}

