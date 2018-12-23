#! @Chapter Calculating centralizer rings

#! @Section Centralizer (commutant) of a representation

#! @Arguments rho

#! @Returns List of standard generators (as a vector space) for the
#! centralizer ring of $\rho(G)$, written in the basis given by <Ref
#! Attr="BlockDiagonalBasis" />.  The matrices are given as a list of
#! blocks.

#! @Description Let $G$ have irreducible representations $\rho_i$ with
#! multiplicities $m_i$. The centralizer has dimension $\sum_i m_i^2$
#! as a $\mathbb{C}$-vector space. This function gives the minimal
#! number of generators required.
DeclareGlobalFunction( "RepresentationCentralizerBlocks" );

#! @Section Useful convenience functions

#! @Arguments rho

#! @Returns List of standard generators (as a vector space) for the
#! centralizer ring of $\rho(G)$.

#! @Description This gives the same result as <Ref
#! Func="RepresentationCentralizerBlocks" />, but with the matrices
#! given in their entirety: not as lists of blocks, but as full
#! matrices.
DeclareGlobalFunction( "RepresentationCentralizer" );
