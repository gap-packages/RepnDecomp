#! @Chapter Computing decompositions of representations

#! @Section Algorithms due to Serre

#! Note: all computation in this section is actually done in the
#! function <Ref Attr="REPN_ComputeUsingSerre" Label="for IsGroupHomomorphism" />,
#!the other functions are wrappers around it.

#! @Arguments rho

#! @Returns List of vector spaces $V_i$, each $G$-invariant and a
#! direct sum of isomorphic irreducibles. That is, for each $i$, $V_i
#! \cong \oplus_j W_i$ (as representations) where $W_i$ is an
#! irreducible $G$-invariant vector space.

#! @Description Computes the canonical decomposition of $V$ into
#! $\oplus_i\;V_i$ using the formulas for projections $V \to V_i$ due
#! to Serre.

#! You can pass in the option `irreps` with a list of irreps of $G$,
#! and this will be used instead of computing a complete list
#! ourselves. If you already know which irreps will appear in $\rho$,
#! for instance, this will save time.
#!
#! <P/>
#! @InsertChunk Example_CanonicalDecomposition
#! <P/>
DeclareGlobalFunction( "CanonicalDecomposition" );

#! @Arguments rho

#! @Returns List of vector spaces $W_j$ such that $V = \oplus_j W_j$
#! and each $W_j$ is an irreducible $G$-invariant vector space.

#! @Description Computes the decomposition of $V$ into irreducible
#! subprepresentations.
#!
#! <P/>
#! @InsertChunk Example_IrreducibleDecomposition
#! <P/>
DeclareGlobalFunction( "IrreducibleDecomposition" );

#! @Arguments rho

#! @Returns List of lists $V_i$ of vector spaces $V_{ij}$ such that $V
#! = \oplus_i \oplus_j V_{ij}$ and $V_{ik} \cong V_{il}$ for all $i$,
#! $k$ and $l$ (as representations).

#! @Description Computes the decomposition of $V$ into irreducible
#! subrepresentations, grouping together the isomorphic
#! subrepresentations.
DeclareGlobalFunction( "IrreducibleDecompositionCollected" );

#! @Arguments rho

#! @Returns A record, in the format described below

#! @Description This function does all of the computation and (since
#! it is an attribute) saves the results. Doing all of the
#! calculations at the same time ensures consistency when it comes to
#! irrep ordering, block ordering and basis ordering. There is no
#! canonical ordering of irreps, so this is crucial.
#!
#! <A>irreps</A> is the complete list of irreps involved in the direct
#! sum decomposition of <A>rho</A>, this can be given in case the
#! default (running Dixon's algorithm) is too expensive, or e.g. you
#! don't want representations over Cyclotomics.
#!
#! The return value of this function is a record with fields:
#!
#! * `basis`: The basis that block diagonalises $\rho$, see <Ref
#!   Func="BlockDiagonalBasisOfRepresentation" />.
#!
#! * `diagonal_rep`: $\rho$, block diagonalised with the basis
#!   above. See <Ref Func="BlockDiagonalRepresentation" />
#!
#! * `decomposition`: The irreducible $G$-invariant subspaces,
#!   collected according to isomorphism, see <Ref
#!   Func="IrreducibleDecompositionCollected" />
#!
#! * `centralizer_basis`: An orthonormal basis for the centralizer
#!   ring of $\rho$, written in block form. See <Ref
#!   Func="CentralizerBlocksOfRepresentation" />
#!
#! Pass the option `parallel` for the computations per-irrep to be
#! done in parallel.
#!
#! Pass the option `irreps` with the complete list of irreps of $\rho$
#! to avoid recomputing this list (could be very expensive)
#!
#! <P/>
#! @InsertChunk Example_REPN_ComputeUsingSerre
#! <P/>
DeclareAttribute( "REPN_ComputeUsingSerre", IsGroupHomomorphism );
