#! @Chapter Centralizer (commutant) rings

#! @Section Finding a basis for the centralizer

#! @Description Let $G$ have irreducible representations $\rho_i$ with
#! multiplicities $m_i$. The centralizer has dimension $\sum_i m_i^2$
#! as a $\mathbb{C}$-vector space. This function gives the minimal
#! number of generators required.
#!
#! <P/>
#! @InsertChunk Example_CentralizerBlocksOfRepresentation
#! <P/>
#! @Arguments rho
#! @Returns List of vector space generators for the centralizer ring
#! of $\rho(G)$, written in the basis given by <Ref
#! Func="BlockDiagonalBasisOfRepresentation" />.  The matrices are
#! given as a list of blocks.
DeclareGlobalFunction( "CentralizerBlocksOfRepresentation" );

#! @Description This gives the same result as <Ref
#! Func="CentralizerBlocksOfRepresentation" />, but with the matrices
#! given in their entirety: not as lists of blocks, but as full
#! matrices.
#!
#! <P/>
#! @InsertChunk Example_CentralizerOfRepresentation
#! <P/>
#! @Returns List of vector space generators for the centralizer ring
#! of $\rho(G)$.
DeclareGlobalFunction( "CentralizerOfRepresentation" );

#! @Section Using the centralizer for computations

#! @Description We require that <A>rho</A> is unitary. Uses the given
#! orthonormal basis (with respect to the inner product $\langle A, B
#! \rangle = \mbox{Trace}(AB^*)$) for the centralizer ring of
#! <A>rho</A> to calculate the sum of the conjugacy class <A>class</A>
#! quickly, i.e. without summing over the class.
#!
#! NOTE: Orthonormality of <A>cent_basis</A> and unitarity of
#! <A>rho</A> are checked. See <Ref Func="ClassSumCentralizerNC" />
#! for a version of this function without checks. The checks are not
#! very expensive, so it is recommended you use the function with
#! checks.
#!
#! <P/>
#! @InsertChunk Example_ClassSumCentralizer
#! <P/>
#! @Arguments rho, class, cent_basis
#! @Returns $\sum_{s \in t^G} \rho(s)$, where $t$ is a representative
#! of the conjugacy class <A>class</A> of $G$.
DeclareGlobalFunction( "ClassSumCentralizer" );

#! @Description The same as <Ref Func="ClassSumCentralizer" />, but
#! does not check the basis for orthonormality or the representation
#! for unitarity.
#!
#! <P/>
#! @InsertChunk Example_ClassSumCentralizerNC
#! <P/>
#! @Arguments rho, class, cent_basis
DeclareGlobalFunction( "ClassSumCentralizerNC" );

DeclareGlobalFunction( "SizesToBlocks" );
