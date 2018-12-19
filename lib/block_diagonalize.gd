#! @Chapter Block diagonalizing representations

#! @Section Finding the correct basis

#! Given a representation $\rho : G \to GL(V)$, it is often desirable
#! to find a basis for $V$ that block diagonalizes each $\rho(g)$ with
#! the block sizes being as small as possible.

#! @Arguments rho

#! @Returns Basis for $V$ that block diagonalizes $\rho$.

#! @Description Let $G$ have irreducible representations $\rho_i$,
#! with dimension $d_i$ and multiplicity $m_i$. The basis returned by
#! this operation gives each $\rho(g)$ as a block diagonal matrix
#! which has $m_i$ blocks of size $d_i \times d_i$ for each $i$.
DeclareAttribute( "BlockDiagonalBasis", IsGroupHomomorphism );

#! @Arguments rho

#! @Returns Representation of $G$ isomorphic to $\rho$ where the
#! images $\rho(g)$ are block diagonalized.

#! @Description This is just a convenience operation that uses <Ref
#! Attr="BlockDiagonalBasis" /> to calculate the basis change matrix
#! and put $\rho$ into a nice form.
DeclareAttribute( "BlockDiagonalRepresentation", IsGroupHomomorphism );

#! @Arguments rho[, irreps]

#! @Returns A record

#! @Description Calculates some values given by other functions in
#! this package, but without summing over G (in case G is very
#! large). Most of the time, G is a symmetric group (or nearly one) so
#! calculation of conjugacy classes is easy. The return value of this
#! function is a record with fields:

#! * basis: same as BlockDiagonalBasis

#! * diagonal_rep: same as BlockDiagonalRepresentation

#! * decomposition: same as IrreducibleDecompositionCollected

#! When I say "the same", I mean up to reordering and isomorphism as
#! representations.
DeclareGlobalFunction( "BlockDiagonalRepresentationFast" );
