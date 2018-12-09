#! @Chapter Computing decompositions of representations

#! @Section Algorithms due to Serre

#! These operations compute various decompositions of a representation
#! $\rho : G \to GL(V)$ where $G$ is finite and $V$ is a
#! finite-dimensional $\mathbb{C}$-vector space.

#! The terms used here are taken from Serre's Linear Representations
#! of Finite Groups.

#! @Arguments rho

#! @Returns List of vector spaces $V_i$, each $G$-invariant and a
#! direct sum of isomorphic irreducibles. That is, for each $i$, $V_i
#! \cong \oplus_j W_i$ (as representations) where $W_i$ is an
#! irreducible $G$-invariant vector space.

#! @Description Computes the canonical decomposition of $V$ into
#! $\oplus_i\;V_i$ using the formulas for projections $V \to V_i$ due
#! to Serre.
DeclareAttribute( "CanonicalDecomposition", IsGroupHomomorphism );

#! @Arguments rho

#! @Returns List of vector spaces $W_j$ such that $V = \oplus_j W_j$
#! and each $W_j$ is an irreducible $G$-invariant vector space.

#! @Description Computes the decomposition of $V$ into irreducible
#! subprepresentations.
DeclareAttribute( "IrreducibleDecomposition", IsGroupHomomorphism );

#! @Arguments rho

#! @Returns List of lists $V_i$ of vector spaces $V_{ij}$ such that $V
#! = \oplus_i \oplus_j V_{ij}$ and $V_{ik} \cong V_{il}$ for all $i$,
#! $k$ and $l$ (as representations).

#! @Description Computes the decomposition of $V$ into irreducible
#! subrepresentations, grouping together the isomorphic
#! subrepresentations.
DeclareAttribute( "IrreducibleDecompositionCollected", IsGroupHomomorphism );
