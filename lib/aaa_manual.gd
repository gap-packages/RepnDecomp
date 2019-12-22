#! @Chapter Introduction

#! @Section Getting started with RepnDecomp

#! This package allows computations of various decompositions of a
#! representation $\rho : G \to GL(V)$ where $G$ is finite and $V$ is
#! a finite-dimensional $\mathbb{C}$-vector space.
#!
#! @Subsection Installation
#!
#! To install this package, refer to the installation instructions in
#! the README file in the source code. It is located here:
#! <URL>https://github.com/gap-packages/RepnDecomp/blob/master/README.md</URL>.

#! @Subsection Note on what is meant by a representation

#! Throughout this documentation, mathematical terminology is used
#! e.g. representation. It is clear what is meant mathematically, but
#! it is not entirely clear what is meant in terms of GAP types - what
#! are you supposed to pass in when I say "pass in a representation".
#! Occasionally I will not even mention what we are passing in and
#! assume the reader knows that <A>rho</A> or $\rho$ refers to a
#! representation.

#! A representation we can use is, in GAP, a homomorphism from a
#! finite group to a matrix group where all matrices have coefficients
#! in a cyclotomic field (`Cyclotomics` is the union of all such
#! fields in GAP). You can check whether something you want to pass is
#! suitable with the function <Ref
#! Attr="IsFiniteGroupLinearRepresentation" Label="for IsGroupHomomorphism"/>.
#!

#! Here's an example of a representation <A>rho</A> in GAP:
#!

#!<Example>
#!gap> G := SymmetricGroup(3);
#!Sym( [ 1 .. 3 ] )
#!gap> images := List(GeneratorsOfGroup(G), g -> PermutationMat(g, 3));
#![ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ],
#!  [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ] ]
#!gap> rho := GroupHomomorphismByImages(G, Group(images));
#![ (1,2,3), (1,2) ] -> [ [ [ 0, 1, 0 ], [ 0, 0, 1 ], [ 1, 0, 0 ] ],
#!  [ [ 0, 1, 0 ], [ 1, 0, 0 ], [ 0, 0, 1 ] ] ]
#!</Example>

#!
#! @Subsection API Overview
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
