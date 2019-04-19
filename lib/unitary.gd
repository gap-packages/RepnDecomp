#! @Chapter Algorithms for unitary representations

#! @Section Unitarising representations

#! @Arguments rho

#! @Returns A record with fields L and unitary_rep such that
#! <A>rho</A> is isomorphic to unitary_rep, differing by a change of
#! basis L.
DeclareGlobalFunction( "UnitaryRepresentation" );

#! @Arguments rho

#! @Returns Whether <A>rho</A> is unitary, i.e. for all $g \in G$,
#! $\rho(g^{-1}) = \rho(g)^*$ (where $^*$ denotes the conjugate
#! transpose).
DeclareGlobalFunction( "IsUnitaryRepresentation" );

#! @Arguments A

#! @Returns a record with two fields, L and D such that $A =
#! L\mbox{diag}(D)L^*$. $D$ is the $1 \times n$ vector which gives the
#! diagonal matrix $\mbox{diag}(D)$ (where <A>A</A> is an $n \times n$
#! matrix).
DeclareGlobalFunction( "LDLDecomposition" );

#! @Section Decomposing unitary representations

#! @Arguments rho

#! @Returns a list of irreps in the decomposition of <A>rho</A>

#! @Description NOTE: this is not implemented yet. Assumes that
#! <A>rho</A> is unitary and uses an algorithm due to Dixon to
#! decompose it into unitary irreps.
DeclareGlobalFunction( "IrreducibleDecompositionDixon" );
