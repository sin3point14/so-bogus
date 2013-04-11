#ifndef BOGUS_NS_NETWON_HPP
#define BOGUS_NS_NETWON_HPP

namespace bogus {

template < typename NSFunction >
class NonSmoothNewton
{
public:
  typedef typename NSFunction::Traits Traits ;
  typedef typename Traits::Scalar Scalar ;
  typedef typename Traits::Vector Vector ;
  typedef typename Traits::Matrix Matrix ;

  NonSmoothNewton( const NSFunction &func, Scalar tol )
  : m_func( func ), m_tol( tol ), m_maxIters( 20 )
  {}

  Scalar solve ( Vector &x ) const ;

private:
  const NSFunction& m_func ;
  Scalar m_tol ;
  unsigned m_maxIters ;

} ;


}

#endif
