---
title: 'RepnDecomp: A GAP package for decomposing linear representations of finite groups'
tags:
  - GAP
  - mathematics
  - groups
  - representations
authors:
  - name: Kaashif Hymabaccus
    orcid: 0000-0001-8581-5804
    affiliation: 1
  - name: Dmitrii Pasechnik
    orcid: 0000-0002-7557-6886
    affiliation: 1
affiliations:
 - name: University of Oxford
   index: 1
date: 31 August 2019
bibliography: paper.bib
---

# Summary

A linear representation of a finite group is a homomorphism from a
finite group G to the group of linear automorphisms of a vector space
V. When studying problems in algebra and combinatorics, it is often
useful to also study associated representations to better understand
the structure of the problem. This is also useful for computation,
since representations allow us to use tools from linear algebra to
solve problems in group theory. A key property of complex
representations of finite groups is that all such representations are
completely reducible, meaning we can decompose them into direct sums
of irreducible representations. When working with matrices, this
corresponds to finding a basis which produces an optimal block
diagonal form of the representation, with the smallest possible
blocks. In cases where there are many small blocks, this can greatly
improve the efficiency of computations done with the matrices.

Currently, while methods for doing these decompositions are known (and
are described in @serre:1977), there are no open-source computer
programs that implement these methods, nor are details on how to
achieve good performance of such an implementation published.

Using the GAP system [@gap:2020], we have produced a package
[RepnDecomp](https://github.com/gap-packages/RepnDecomp) that provides
the following functions:

* Decompose a representation into a direct sum of irreducible
  representations
* Compute the associated basis that gives an optimal block
  diagonalisation
* Determine whether two representations are isomorphic and compute the
  isomorphism
* Compute a unitary representation isomorphic to a given
  representation
* Compute the centraliser of the representation (the vector space of
  matrices that commute with the matrices of the representation)

Our package deals exclusively with the case where V is a
finite-dimensional $\mathbb{C}$-vector space with linear automorphisms
represented as square matrices. In fact, we only consider cases where
the matrix coefficients are cyclotomic numbers - complex numbers in
the $\mathbb{Q}$-vector space spanned by all powers of all roots of
unity. Our methods are not specific to cyclotomic numbers, but the GAP
system only has facilities to compute with cyclotomics and not general
complex numbers.

This package has been applied to improve the block diagonalisation of
matrices involved in a semidefinite program for computing bounds on
the crossing number of complete graphs [@deklerk:2007]. Our package
can also be applied to other problems mentioned by @deklerk:2007,
including the computation of bounds for the Lovász $\vartheta$ (and
related $\vartheta'$) numbers for graphs, and the truss topology
design problem described in @kanno:1970 for trusses with suitable
symmetries.

The algorithms used to implement this package are not all
original. One algorithm for decomposing representations is original
work. The other is based on formulas described in @serre:1977 ---
however, we implement apparently novel speedups to make it feasible
for groups of large order. Our method for unitarising a representation
is based on @dixon:1970, again, with apparently novel speedups.

The development of the package was part of Kaashif Hymabaccus'
Master's thesis [@hymabaccus:2019].

# Acknowledgements

Kaashif Hymabaccus thanks the Collaborative Computational Project
(CCP) in the area of Computational Discrete Mathematics, supported by
the EPSRC (EP/M022641/1) for supporting his participation in a GAP
workshop.

Dmitrii Pasechnik was partially funded by the OpenDreamKit Horizon
2020 European Research Infrastructures project (#676541).

# References
