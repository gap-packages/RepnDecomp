# Computing decompositions of representations

## Project Description

This project aims to investigate and record algorithmic details and to
produce a program (using GAP and/or Sagemath) for computing the
decomposition of a representation ρ, of a finite group G over the
complex numbers into irreducibles, as well as the corresponding
decomposition of the centraliser of R.  Currently, while methods for
doing this are known (cf.  Serre's book "Linear Representations of
Finite Groups"), there are no open-source computer programs that
implement these methods, nor are details on how to achieve good
performance of such an implementation published.

This program will be useful in, for example, semidefinite programming
and optimisation/feasibility problems involving coding theory, graph
theory, algebraic geometry, combinatorics and more (see
https://arxiv.org/abs/1007.2905 for more examples of possible
applications). Specifically, it allows to achieve substantial
reductions in the dimension of these problems; potentially known
results, e.g. on upper bounds on sizes of nonlinear codes, could be
improved with the help of the program.

## Testing

To run the tests, make sure you have all needed packages installed (GRAPE).

    $ gap tst/testall.g

This will run all tests and (hopefully) pass.

## Documentation

There's a GAPDoc documentation book hosted
[here](https://kaashif.gitlab.io/decomp/chap0.html). This is generated
from the source files and comments in the `lib/` directory of this
repo, so you can also look there for the same information.

There are also some examples in the `examples` directory.

## Paper

This package has been submitted as a paper to the
[Journal of Open Source Software](https://joss.theoj.org/), the paper
can be found at `paper.md`.
