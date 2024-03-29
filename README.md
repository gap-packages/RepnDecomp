[![Build Status](https://github.com/gap-packages/RepnDecomp/workflows/CI/badge.svg?branch=master)](https://github.com/gap-packages/RepnDecomp/actions?query=workflow%3ACI+branch%3Amaster)
[![Code Coverage](https://codecov.io/github/gap-packages/RepnDecomp/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/RepnDecomp)

# Computing decompositions of representations

## Overview

The main function of this package is to compute the decomposition of a
representation ρ, of a finite group G over the complex numbers into
irreducibles, as well as the corresponding decomposition of the
centraliser of R.

While methods for doing this were well known before this package was
written (cf.  Serre's book "Linear Representations of Finite Groups"),
there were no open-source computer programs that implemented these
methods, nor were details on how to achieve good performance of such
an implementation published.

This package is useful in, for example, semidefinite programming and
optimisation/feasibility problems involving coding theory, graph
theory, algebraic geometry, combinatorics and more (see
https://arxiv.org/abs/1007.2905 for more examples of possible
applications). Specifically, it allows to achieve substantial
reductions in the dimension of these problems; potentially known
results, e.g. on upper bounds on sizes of nonlinear codes, could be
improved with the functions implemented in this package.

This package was written as part of my (Kaashif Hymabaccus's) Master's
degree at the University of Oxford, supervised by Dmitrii Pasechnik.

## Installation

Make sure that, when you install GAP, you installed all of the
packages! Our package uses GRAPE and IO for some functions.

### Latest version included with GAP

If you have version 4.11.0 or later of the
[GAP system](https://www.gap-system.org/Download/index.html)
installed, you do not have to install RepnDecomp manually since it is
already distributed with GAP.

### Latest released version

If you would like the latest released version of RepnDecomp, and there
has not yet been a release of GAP including it, then you can download
the latest release
[here](https://gap-packages.github.io/RepnDecomp/).

Create the directory `~/.gap/pkg`, which will contain your local
packages and extract the archive you downloaded into it. For example:

    $ mkdir -p ~/.gap/pkg
    $ tar -C ~/.gap/pkg -xvzf RepnDecomp-1.1.0.tar.gz

### Development version

If you would like to install the latest code directly from git master
(unsuitable for anything but development of this package), then you
can clone this repo directly into `~/.gap/pkg`:

    $ mkdir -p ~/.gap/pkg
    $ cd ~/.gap/pkg
    $ git clone https://github.com/gap-packages/RepnDecomp.git

### Post-installation

Now, you can run GAP however you like, load the package and use the
functions provided:

```
$ gap
<some output>
gap> LoadPackage("RepnDecomp");
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

Make sure your current directory is this repo, for example:

    $ git clone https://github.com/gap-packages/RepnDecomp.git
    $ cd RepnDecomp

First, generate the tests from the documentation:

    $ gap -q < makedoc.g

To run the tests, make sure you have all needed packages installed
(GRAPE and also IO if you want to compute in parallel).

    $ gap tst/testall.g

This will run all tests and (hopefully) pass. In order for the tests
to be as useful as possible to me i.e. catch as many bugs as possible,
there is a lot of randomness in them. This means that sometimes the
tests pick a pathologically bad example to decompose which causes the
tests to hang forever. The tests will be fixed to be more
deterministic.

## Documentation

There's a GAPDoc documentation book hosted
[here](https://gap-packages.github.io/RepnDecomp/doc/chap0.html). This
is generated from the source files and comments in the `lib/`
directory of this repo, so you can also look there for the same
information.

There are also some examples in the `examples` directory, which are
embedded in the manual, but the most complete examples are in the
`old_tst` directory. These are the old, poorly documented, but
essentially complete tests. They are being converted into nicer
looking tests that can be embedded into the manual, but this could
take some time.

## Paper

This package is described in a freely available [paper](https://joss.theoj.org/papers/10.21105/joss.01835) in
[Journal of Open Source Software](https://joss.theoj.org/).

## Contributing

Open a pull request or issue in this repository. You can also email me
directly, I'll be more likely to notice an email.
