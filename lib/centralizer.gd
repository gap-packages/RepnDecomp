#! @Arguments rho
#! @Returns list of standard generators for centralizer of rho(G) in block form
#! @Description
#!   Decomposes rho into irreducibles, block diagonalizes rho, then
#!   uses the block sizes to calculate the centralizer C of rho. The
#!   generators of C are decomposed into blocks, each block
#!   corresponding to an isomorphism class of irrep in the decomposition
#!   of rho.
DeclareGlobalFunction( "RepresentationCentralizerBlocks" );

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
#! @Returns list of standard generators for centralizer of rho(G)
#! @Description
#!   Decomposes character of rho into sum of irreducible characters,
#!   uses that to calculate basis for centralizer of rho.
DeclareGlobalFunction( "RepresentationCentralizer" );
