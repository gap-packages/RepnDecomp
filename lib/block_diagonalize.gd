#! @Arguments rho
#! @Returns isomorphic representation where images are block diagonal
#! @Description
#!   Calculates the matrix A such that for any g in G, A^-1 * rho(g) *
#!   A is block diagonal with each block corresponding to an
#!   irreducible in the decomposition of rho. The blocks are ordered
#!   so the isomorphic irreps' blocks are adjacent. Returns a
#!   representation isomorphic to rho that uses this basis change
#!   matrix.
DeclareAttribute( "BlockDiagonalRepresentation", IsGroupHomomorphism );

DeclareAttribute( "BlockDiagonalBasis", IsGroupHomomorphism );
