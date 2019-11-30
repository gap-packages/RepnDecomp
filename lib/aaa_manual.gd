#! @Chapter Introduction

#! @Section Getting started with RepnDecomp

#! This package allows computations of various decompositions of a
#! representation $\rho : G \to GL(V)$ where $G$ is finite and $V$ is
#! a finite-dimensional $\mathbb{C}$-vector space.
#!
#! To install this package, refer to the installation instructions in
#! the README file in the source code. It is located here:
#! <URL>https://gitlab.com/kaashif/decomp/blob/master/README.md</URL>.
#!
#! The algorithms implemented can be divided into two groups: methods
#! due to Serre from his book Linear Representations of Finite Groups,
#! and original methods due to the authors of this package.
#!
#! The default is to use the algorithms due to Serre. If you pass the
#! option `method := "alternate"` to a function, it will use the
#! alternate method. Passing the option `parallel` will try to compute
#! in parallel as much as possible. See the individual functions for
#! options you can pass.
#!
#! The main functions implemented in this package are:
#!
#! For decomposing representations into canonical and irreducible
#! direct summands:
#!
#! * <Ref Func="CanonicalDecomposition" />
#! * <Ref Func="IrreducibleDecomposition" />
#! * <Ref Func="IrreducibleDecompositionCollected" />
#!
#! For block diagonalising representations:
#!
#! * <Ref Func="BlockDiagonalBasisOfRepresentation" />
#! * <Ref Func="BlockDiagonalRepresentation" />
#!
#! For computing centraliser rings:
#!
#! * <Ref Func="CentralizerBlocksOfRepresentation" />
#! * <Ref Func="CentralizerOfRepresentation" />
#!
#! For testing isomorphism and computing isomorphisms (intertwining
#! operators) between representations:
#!
#! * <Ref Func="LinearRepresentationIsomorphism" />
#! * <Ref Func="AreRepsIsomorphic" />
#! * <Ref Func="IsLinearRepresentationIsomorphism" />
#!
#! For testing unitarity of representations and the unitarisation of
#! representations:
#!
#! * <Ref Func="UnitaryRepresentation" />
#! * <Ref Func="IsUnitaryRepresentation" />

#! @Chapter Decomposing representations into irreducibles

#! @Chapter Centraliser (commutant) rings

#! @Chapter Isomorphisms between representations

#! @Chapter Algorithms for unitary representations

#! @Chapter Miscellaneous useful functions
