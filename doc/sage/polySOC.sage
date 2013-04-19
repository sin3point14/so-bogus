var('u, r, t, CT, ST')
var ( 'w00, w01, w02, w11, w12, w22, b0, b1, b2' ) 
var('mu')

eq1 = u == ( w00 + w01 * mu * CT + w02 * mu * ST )*r + b0

eq2 = -u * CT / mu  - ( w01 + w11 * mu * CT + w12 * mu * ST )*r ==  b1
eq3 = -u * ST / mu  - ( w02 + w12 * mu * CT + w22 * mu * ST )*r ==  b2

eq2bis = eq2.substitute(u = eq1.right()).expand()
eq3bis = eq3.substitute(u = eq1.right()).expand()

eq2c = eq2bis.add_to_both_sides( b0*CT/mu )
eq3c = eq3bis.add_to_both_sides( b0*ST/mu )

eq = eq2c * eq3c.right() - eq3c * eq2c.right()
assume( r>0 )
assume( mu>0 )

eqc = ( eq * mu / r).expand()

var ( 'D' )
var ( 't' )

eqt = eqc.substitute( ST = 2*t/D ).substitute( CT = (1 - t*t)/D ) * D * D
eqt = eqt.expand().substitute( D = 1+t*t ).expand()

P = eqt.left().poly(t)
Pcoeffs = P.coeffs(t)

for i in range(0,4):
  print "coeffs[",i,"]= ", Pcoeffs[i][0].poly(mu), ";" 

print solve(eq2bis, r )

print "--- Numerical test ---"

d_mu = 1
d_b0 = -1 
d_b1 = 1 
d_b2 = 1
d_w00 = 1
d_w11 = 2
d_w22 = 2
d_w01 = 0
d_w02 = 0
d_w12 = 0

d_mu = 1
d_b0 = -1 
d_b1 = 1 
d_b2 = 1
d_w00 = 1
d_w11 = 2
d_w22 = 2
d_w01 = 0
d_w02 = 0
d_w12 = 0

data = { mu:d_mu, b0:d_b0, b1:d_b1, b2:d_b2, w00:d_w00, w11:d_w11, w22:d_w22, w01:d_w01, w02:d_w02, w12:d_w12 } 
Ptest = P.substitute( data )
Peq = Ptest == 0

t0 = find_root( Peq, -1.e6, 1.e6 )
print "tan(theta/2)= ", t0

theta0 = 2*arctan( t0 )
CT0=cos(theta0)
ST0=sin(theta0)

eq2test = eq2bis.substitute( data ).substitute( {CT:CT0, ST:ST0 } )
r0 = find_root( eq2test, -1.e6, 1.e6 )
#r0 = eq2test.roots()[0][0]
eq1test = eq1.substitute( data ).substitute( {CT:CT0, ST:ST0, r:r0 } )
u0 = find_root( eq1test, -1.e6, 1.e6 )
#u0 = eq1test.roots()[0][0]


r1 = d_mu * r0 * CT0
r2 = d_mu * r0 * ST0

u1 = -u0 * CT0 / d_mu
u2 = -u0 * ST0 / d_mu

u0h = d_b0 + d_w00 * r0 + d_w01 * r1 + d_w02 * r2
u1h = d_b1 + d_w01 * r0 + d_w11 * r1 + d_w12 * r2
u2h = d_b2 + d_w02 * r0 + d_w12 * r1 + d_w22 * r2

d0 = (u0h - u0)
d1 = (u1h - u1)
d2 = (u2h - u2)

print "rN= ", r0, " uN= ", u0
print "|rT|= ", sqrt(r1*r1 + r2*r2), " |uT|= ", sqrt(u1*u1 + u2*u2)
print "r.dot(u) = ", u0*r0 + u1*r1 + u2*r2
print "err = ", d0, d1, d2

