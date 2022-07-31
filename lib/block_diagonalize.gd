#! @Chapter Computing decompositions of representations

#! @Section Block diagonalizing

#! Given a representation $\rho : G \to GL(V)$, it is often desirable
#! to find a basis for $V$ that block diagonalizes each $\rho(g)$ with
#! the block sizes being as small as possible. This speeds up matrix
#! algebra operations, since they can now be done block-wise.

#! @Arguments rho

#! @Returns Basis for $V$ that block diagonalizes $\rho$.

#! @Description Let $G$ have irreducible representations $\rho_i$,
#! with dimension $d_i$ and multiplicity $m_i$. The basis returned by
#! this operation gives each $\rho(g)$ as a block diagonal matrix
#! which has $m_i$ blocks of size $d_i \times d_i$ for each $i$.
DeclareGlobalFunction( "BlockDiagonalBasisOfRepresentation" );

#! @Arguments rho

#! @Returns Representation of $G$ isomorphic to $\rho$ where the
#! images $\rho(g)$ are block diagonalized.

#! @Description This is just a convenience operation that uses <Ref
#! Attr="BlockDiagonalBasisOfRepresentation" /> to calculate the basis
#! change matrix and applies it to put $\rho$ into the block
#! diagonalised form.
DeclareGlobalFunction( "BlockDiagonalRepresentation" );

DeclareGlobalFunction( "BasisChangeMatrixSimilar" );
