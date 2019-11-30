# Computing decompositions of representations

## Project Description

This project aims to investigate and record algorithmic details and to
produce a program (using GAP and/or Sagemath) for computing the
decomposition of a representation ρ, of a finite group G over the
complex numbers into irreducibles, as well as the corresponding
decomposition of the centraliser of R. Currently, while methods for
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

## Installation

First, install GAP following the instructions
[here](https://www.gap-system.org/Download/index.html). Then, create a
directory `~/.gap/pkg`, which will contain your local packages and
clone this into it. Commands to run:

    $ mkdir -p ~/.gap/pkg
    $ cd ~/.gap/pkg
    $ git clone https://gitlab.com/kaashif/decomp.git RepnDecomp

Make sure that, when you installed GAP, you installed all of the
packages! Our package uses GRAPE and IO for some functions.

Now, you can run GAP however you like, load the package and use the
functions provided:

```
$ gap
<some output>
gap> gap> LoadPackage("RepnDecomp");
───────────────────────────────────────────────────────────────────────────────
Loading  GRAPE 4.8.2 (GRaph Algorithms using PErmutation groups)
by Leonard H. Soicher (http://www.maths.qmul.ac.uk/~lsoicher/).
Homepage: https://gap-packages.github.io/grape
Report issues at https://github.com/gap-packages/grape/issues
───────────────────────────────────────────────────────────────────────────────
───────────────────────────────────────────────────────────────────────────────
Loading  RepnDecomp 0.1 (Decompose representations of finite groups into irreducibles)
by Kaashif Hymabaccus (https://kaashif.co.uk).
with contributions by:
   Dmitrii Pasechnik.
Homepage: http://gitlab.com/kaashif/decomp/
───────────────────────────────────────────────────────────────────────────────
true
gap> A := IdentityMat(5);
[ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, 0, 0, 1, 0 ],
  [ 0, 0, 0, 0, 1 ] ]
gap> B := LDLDecomposition(A);
rec( D := [ 1, 1, 1, 1, 1 ],
  L := [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
      [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] )
```

Where `LDLDecomposition` is a function provided by this package.

## Testing

To run the tests, make sure you have all needed packages installed
(GRAPE and also IO if you want to compute in parallel).

    $ gap tst/testall.g

This will run all tests and (hopefully) pass.

## Documentation

There's a GAPDoc documentation book hosted
[here](https://kaashif.gitlab.io/decomp/chap0.html). This is generated
from the source files and comments in the `lib/` directory of this
repo, so you can also look there for the same information.

There are also some examples in the `examples` directory, but the most
complete set of examples are in the `tst` directory - the unit tests
themselves.

## Paper

This package has been submitted as a paper to the
[Journal of Open Source Software](https://joss.theoj.org/), the paper
can be found at `paper.md`.
