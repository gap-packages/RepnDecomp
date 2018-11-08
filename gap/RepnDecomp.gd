#
# RepnDecomp: Decompose representations of finite groups into irreducibles
#
# Declarations
#

#! @Arguments rho
#! @Returns a list of vector spaces that direct sum to the given representation
#! @Description
#!   Decomposes <A>rho</A> into a direct sum of subrepresentations,
#!   each irreducible.
DeclareGlobalFunction( "DecomposeRepresentationIrreducible" );

#! @Arguments rho
#! @Returns a list of vector spaces that direct sum to the given representation
#! @Description
#!   Decomposes <A>rho</A> into the canonical representation. Can be
#!   seen as (although not implemented this way) doing the full
#!   decomposition into irreducibles then collecting the isomorphic
#!   representations.
DeclareGlobalFunction( "DecomposeRepresentationCanonical" );

#! @Arguments rho
#! @Returns isomorphic representation where images are block diagonal
#! @Description
#!   Calculates the matrix A such that for any g in G, A^-1 * rho(g) *
#!   A is block diagonal with each block corresponding to an
#!   irreducible in the decomposition of rho. The blocks are ordered
#!   so the isomorphic irreps' blocks are adjacent. Returns a
#!   representation isomorphic to rho that uses this base change
#!   matrix.
DeclareGlobalFunction( "BlockDiagonalizeRepresentation" );

#! @Arguments rho
#! @Returns list of standard generators for centralizer of rho(G)
#! @Description
#!   Decomposes character of rho into sum of irreducible characters,
#!   uses that to calculate basis for centralizer of rho.
DeclareGlobalFunction( "RepresentationCentralizer" );

#! @Arguments rho
#! @Returns list of decomposed standard generators for centralizer of
#!   rho(G) in block form
#! @Description
#!   Decomposes rho into irreducibles, block diagonalizes rho, then
#!   uses the block sizes to calculate the centralizer of rho. This is
#!   the same as RepresentationCentralizerBlocks, but with all
#!   identity blocks compressed into 1x1 cells.
#!
#!   This shows that C is isomorphic to C' a direct sum of m(X_k) x
#!   m(X_k) matrices.
DeclareGlobalFunction( "RepresentationCentralizerDecomposed" );

#! @Arguments rho
#! @Returns list of standard generators for centralizer of rho(G) in block form
#! @Description
#!   Decomposes rho into irreducibles, block diagonalizes rho, then
#!   uses the block sizes to calculate the centralizer C of rho. The
#!   generators of C are decomposed into blocks, each block
#!   corresponding to an isomorphism class of irrep in the decomposition
#!   of rho.
DeclareGlobalFunction( "RepresentationCentralizerBlocks" );
