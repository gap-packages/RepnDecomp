#! @Chapter Paralellised functions

#! @Section Decompositions

#! These functions are copies of other functions in the package, but
#! made to run in parallel. Usually, the trick is to rearrange the
#! formulas so that we can run independent calculations per irrep in
#! the list of (relevant) irreps of our group $G$.

#! @Arguments rho, num_jobs[, irreps]

#! @Returns List of lists $V_i$ of vector spaces $V_{ij}$ such that $V
#! = \oplus_i \oplus_j V_{ij}$ and $V_{ik} \cong V_{il}$ for all $i$,
#! $k$ and $l$ (as representations).

#! @Description This function uses the same algorithm as <Ref
#! Func="IrreducibleDecompositionCollectedHybrid" />, but runs the
#! per-irrep calculation to break down the canonical summands in
#! parallel. Uses <A>irreps</A> as the list of relevant irreps, if
#! given.
DeclareGlobalFunction( "IrreducibleDecompositionCollectedParallel" );

#! @Arguments rho, num_jobs[, irreps]

#! @Returns The same result as <Ref
#! Func="BlockDiagonalRepresentationFast" />.
DeclareGlobalFunction( "BlockDiagonalRepresentationParallel" );
