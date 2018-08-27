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
