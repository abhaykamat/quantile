quantile
========

quantile.h contains C++ code for an arbitrary precision
evaluation of the quantile function of the normal distribution.

test.cpp contains test/example code

Code is geared towards precision, and not speed. Typical use
is to verify how well the fast machine-precision approximations
of the quantile work as advertized.

It's unclear how well it actually behaves deep in the tails
(i.e.  when variable is extremely close to zero or one).

mpreal.h is stolen from here :

    http://www.holoborodko.com/pavel/mpfr/

Code is in the public domain.

Package requires mpfr

