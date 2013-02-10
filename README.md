<!-- this file uses MarkDown syntax,
    see: http://daringfireball.net/projects/markdown/ -->

This repository hosts sample code implementing the
[Raghavendra][1]-[Fliege][3] algorithm for solving linear systems.


Available code
--------------

The code is grouped into directories according to the implementation
language.  All example programs run two tests in sequence: solve a
linear system `Ax = b` with `A` the identiy matrix, and solve a linear
system where `A` is a random matrix with coefficients drawn uniformly
from an interval.


References
----------

* The original [draft paper by Raghavendra][1] describing the new
  solver for linear systems on finite fields:
  <http://www.eecs.berkeley.edu/~prasad/linsystems.pdf>

* [R. J. Lipton's blog post][2] discussing Raghavendra's algorithm:
  <http://rjlipton.wordpress.com/2012/08/09/a-new-way-to-solve-linear-equations/>

* The [arXiv paper by J. Fliege][3] providing a modified version of
  the algorithm that works over the reals:
  <http://arxiv.org/pdf/1209.3995v1>

[1]: http://www.eecs.berkeley.edu/~prasad/linsystems.pdf
[2]: http://rjlipton.wordpress.com/2012/08/09/a-new-way-to-solve-linear-equations/
[3]: http://arxiv.org/pdf/1209.3995v1
